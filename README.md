# Paired FASTQ Service（配对与质控重校）

测序湖入湖工具，包含两个阶段：

1. **配对**：接收顺序被打乱的双端未压缩 FASTQ（R1/R2），按读段号重新
   配对，按 R1 顺序交付两份文件，原文逐字节保留。
2. **对照质量分重校**：对配对成功的任务，按对照与参考的错配率改写双端
   FASTQ 的质量行，打包为 ZIP 供变异分析。

FastAPI + SQLite 实现，单文件上限 256 MiB，进程内存上限 256 MiB，配对
流程可处理总量超过内存的输入。

## 配对规则

- 读段号取头行首个空白前的字段，去掉开头 `@`，再去掉末尾的 `/1` 或 `/2`
  （例如 `@read-7/1 sample=x` → `read-7`）。
- 输出按 **R1 的原始顺序**排列；R2 按相同读段号对齐到同一位置。
- 四行记录原文（含头行描述、换行风格）原样写入结果，不做任何改写。

## 失败即整次失败，不开放部分结果

任一情况发生，任务整体标记为 `failed`/`cancelled`，已写入的暂存全部删除，
下载接口返回 409：

- 同一文件内读段号重复（含 `/1`、`/2` 归一化后重复）
- 孤端（某读段号只出现在一侧）
- 格式错误：非四行式、头行不以 `@` 开头、第三行不以 `+` 开头、空文件
- 序列与质量字符串长度不等

## 对照质量分重校（`POST /tasks/{id}/recalibrate`）

对一个**已成功**的配对任务发起同步请求，服务读取该任务发布的 R1/R2 结果，
按下列规则统计错配并改写质量行，返回内含 `R1.fastq`、`R2.fastq` 的 ZIP。

### 输入约束

- `reference`：至多 10 000 位的 ACGT 参考序列（非空、不含其他字符）。
- `manifest`：每个读对一条，字段为 `read_id`、`r1_start`/`r2_start`
  （1 基起点，整数 ≥ 1）、`r1_strand`/`r2_strand`（`+`/`-`）。
- 读对总数 ≤ 2000；每端读长 40～200 位且只含 ACGT。
- 质量为 Phred+33，分值限定 0～40（字节 `!`～`I`）；无插缺（不解析 CIGAR，
  读段按连续坐标覆盖参考）。
- 反向链（`-`）：统计前先对序列做**反向互补**，并将质量串**倒转**。
- 起点约定为该读覆盖的最左（最小）参考坐标；无论链向，窗口均为
  `[start, start + 读长 - 1]`，越界即整次报错。

### 统计与评分

- 先按**原质量分**统计与参考不符的碱基。
- 同一读对的两端覆盖同一坐标时：
  - 碱基相同 → 该坐标只计一次，取**质量较高者**的碱基参与比对；
  - 碱基不同 → 该坐标两端都不计（冲突位点）。
- 有效位点（去重后实际参与统计的坐标数）达到 **20 个**时，按

  ```
  Q = round_half_up(-10 * log10((错配数 + 1) / (有效位点数 + 2)))
  ```

  计算并四舍五入（半数向上），再裁剪到 0～40；不足 20 个有效位点时
  **沿用原质量分**，不改写。
- 改写时仅替换质量行的字符，头行、`+` 行、序列行及换行风格（含 CRLF）
  等其余原文字节保持不变；R1 顺序不变，R2 按读段号对齐到同一顺序。

### 失败即整次失败

以下任一情况返回错误响应且**不产生任何 ZIP**（更没有部分结果），原任务
文件保持不变：

- 任务号不存在（404）或任务未成功 / 结果缺件（409）；
- manifest 漏号、读段号重复或出现 FASTQ 中不存在的未知读号；
- 起点越界、链向取值非法、参考超长或含非 ACGT 字符；
- FASTQ 后段损坏（非四行式、头行/`+` 行非法、截断、尾部多余字节等）、
  读长越界、质量超出 0～40、读对数超过 2000。

响应头同时回传统计量，便于核对：`X-Effective-Sites`、`X-Mismatches`、
`X-Recalibrated-Score`（沿用原分时为空）、`X-Score-Applied`（`0`/`1`）。

### 内存

重校输入有界（2000 对 × 200 bp，两份 FASTQ 各约数百 KiB，参考 ≤ 10 KB），
记录以小对象保存在内存中，ZIP 全量在 `BytesIO` 中组装；实测 2000 对 ×
200 bp 规模峰值 RSS 约 60 MiB，远低于 256 MiB 上限。

## 内存与大输入

- 上传经 `SpooledTemporaryFile` 落盘（超过 1 MiB 即写临时文件），分块复制，
  单文件流式校验 256 MiB 上限，超限返回 413。
- 解析按 512 条/批写入 SQLite（WAL、小 page cache、临时表走文件
  `temp_store=FILE`、关闭 mmap），任意时刻内存里只有一个批次。
- 配对发布使用 `JOIN ... ORDER BY a.ord` 流式游标，分批 `fetchmany` 写出。
- 因此输入总量可以远超 256 MiB；自测以约 393 MiB 总量（两份各 196.5 MiB、
  各 200 万条）验证，峰值 RSS 约 19 MiB（见下文自测结果）。

## 任务与取消语义

- `POST /tasks` 上传后立即返回 `202`，后台守护线程处理。
- `GET /tasks/{id}` 查询状态、阶段（upload/parse/validate/publish/done）与
  0–1 进度。
- `POST /tasks/{id}/cancel`：置取消标记后**立即返回**，不等待配对结束；
  工作线程在下一个批次边界中止。取消与完成之间用条件 UPDATE 仲裁：
  - 取消先到 → 任务最终 `cancelled`，结果永不发布，暂存目录删除；
  - 完成先到 → 取消请求得到 `409`，明确告知任务已完成、不可取消，结果保留。
- 取消后可立刻创建新任务；任务目录以 UUID 隔离，清理只删除本任务目录，
  不会误删新任务文件。
- 服务重启（启动钩子）：所有非终态任务标记 `failed`（原因
  `interrupted by service restart`）并清理暂存与记录；`succeeded` 任务的
  结果文件保留，仍可下载。

## API

| 方法 | 路径 | 说明 |
| --- | --- | --- |
| POST | `/tasks` | multipart 字段 `r1`、`r2`，返回任务对象（202） |
| GET  | `/tasks/{id}` | 状态/进度/消息/下载链接 |
| POST | `/tasks/{id}/cancel` | 请求取消（立即返回） |
| GET  | `/tasks/{id}/download/r1` | 下载配对后的 R1 |
| GET  | `/tasks/{id}/download/r2` | 下载配对后的 R2 |
| POST | `/tasks/{id}/recalibrate` | 对成功任务做质控重校，返回 ZIP |
| GET  | `/health` | 健康检查 |

启动后打开交互式文档：<http://localhost:8000/docs>，可直接在页面上演示
上传、取消、下载全流程。

## 运行（Docker Compose，内存上限 256 MiB）

```bash
docker compose up --build
# 打开 http://localhost:8000/docs
```

`docker-compose.yml` 中 `mem_limit: 256m` 即交付要求的内存上限；数据保存在
命名卷 `fastq-data`（容器内 `/data`，含 SQLite 与各任务目录），重启后成功
结果仍可下载。

## 本地直接运行（开发/自测）

```bash
python3 -m venv --without-pip .venv
.venv/bin/python tools/get-pip.py
.venv/bin/pip install -r requirements.txt
FASTQ_DATA_DIR=./data .venv/bin/uvicorn app.main:app --port 8000
```

## 演示脚本

```bash
./scripts/demo.sh                 # 需要一个运行在 localhost:8000 的服务
BASE_URL=http://localhost:8000 ./scripts/demo.sh
```

脚本会生成一对乱序 FASTQ，演示上传 → 轮询进度 → 下载校验，并额外演示一次
取消与一次错误输入。

## 自测

```bash
# 短测试（功能/取消/重启语义 + 取消与完成的仲裁 + 质控重校）
bash scripts/run_tests.sh
# 或分别运行
.venv/bin/python -I tests/test_api.py            # 30 项：配对、整次失败、取消、重启恢复
.venv/bin/python -I tests/test_recalibrate.py    # 41 项：重校算法、改写、整次报错
.venv/bin/python -I tests/test_cancel_race.py    # 9 项：取消/完成条件 UPDATE 仲裁

# 内存压力测试：两份各 196.5 MiB（合计约 393 MiB，超过 256 MiB 内存上限），
# 子进程实测峰值 RSS 约 19 MiB，并逐条校验配对顺序与原文
.venv/bin/python -I tests/test_memory.py 2000000
```

## 目录结构

```
app/
  config.py       # 上限、路径、批大小等配置（环境变量可覆盖）
  db.py           # SQLite schema 与任务/记录访问
  fastq.py        # 流式解析、校验、JOIN 配对发布
  recalibrate.py  # 重校：FASTQ 校验、错配统计、评分、质量行改写、ZIP
  worker.py       # 后台线程、取消仲裁、启动恢复、任务目录清理
  main.py         # FastAPI 路由、有界流式上传、重校端点
tests/            # 功能自测与内存压力测试
scripts/demo.sh
Dockerfile, docker-compose.yml, requirements.txt
```

## 遗留项 / 边界说明

- 仅接受未压缩 FASTQ；`.gz` 等压缩格式不在范围内。
- 服务为单实例 + 单 SQLite 文件设计；多实例横向扩展需要换用共享存储/队列。
- 没有内置任务 TTL 清理，成功结果会一直保留在数据卷中，需外部运维定期清理。
- 开发沙箱内未安装 Docker，`docker compose up` 未实跑；compose/Dockerfile 仅做
  了 YAML 与结构校验，应用本身已用等价的本地 uvicorn 命令完整验证（含
  `/docs`、上传、取消、下载、重启恢复）。在装有 Docker 的环境执行
  `docker compose up --build` 即可，`mem_limit: 256m` 即内存上限。
