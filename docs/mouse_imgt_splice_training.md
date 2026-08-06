# 小鼠 OpenSpliceAI：NCBI + IMGT 重新训练中文指南

## 1. 这套代码解决什么问题

作者原始流程主要从 NCBI GFF3 的普通转录本和相邻 exon 推导剪接标签：

- `0`：不是剪接位点
- `1`：acceptor（受体）
- `2`：donor（供体）

但是小鼠免疫球蛋白和 TCR 基因座有特殊情况。例如：

- 单外显子的 `IGHJ` 片段，作者的“相邻 exon”算法不能产生其末端 donor；
- `IGHC` 第一个恒定区 exon 的 acceptor，上游 exon 可能属于重排后的 J 区，不能只靠同一个 `C_gene_segment` 内部 exon 推导；
- NCBI GFF3 可能已经明确标注一部分免疫剪接位点，但这些 `J/C_gene_segment` 不一定会进入作者的 `protein-coding` 基线 HDF5；
- IMGT/LIGM-DB 会明确标注 `DONOR-SPLICE`、`ACCEPTOR-SPLICE`、`INT-DONOR-SPLICE` 和 `INT-ACCEPTOR-SPLICE`。

本分支采用以下优先级：

1. 保留作者原始 NCBI 普通转录本基线；
2. NCBI 已明确注释、但没有真正进入基线 HDF5 的 IG/TR 位点，补入训练集；
3. NCBI 没覆盖的 IMGT functional 位点，再补入训练集；
4. 原始 NCBI test 集完全不变，避免测试集污染。

> OpenSpliceAI 学的是每个碱基属于 donor、acceptor 还是非剪接位点。  
> 它不会直接学习一个“IGHJ1 与 IGHM 配对”的类别。正确做法是让模型看到明确的 J donor 和 C acceptor，以及它们周围的真实序列上下文。

---

## 2. 当前云端数据

数据已经发布在你的 fork 的 GitHub Release：

```text
仓库：Ethan-Fu-2000/OpenSpliceAI
Release：mouse-immune-data-20260806-31072225004
组装版本：GCF_000001635.27（GRCm39）
```

Release 包含：

- NCBI GRCm39 genome FASTA
- NCBI GRCm39 GFF3
- NCBI GRCm39 GBFF
- NCBI sequence report / assembly report
- IMGT/GENE-DB release 信息
- IMGT/LIGM-DB `imgt.dat.Z`
- `manifest.json`
- 每个文件和压缩包的 SHA-256

大型数据没有进入 Git 历史，而是作为 Release 附件保存。

---

## 3. 安装环境

建议单独创建环境，不要污染 base：

```bash
umask 0002

mamba create -n openspliceai_mouse \
  -c conda-forge -c bioconda \
  python=3.11 \
  git \
  gh \
  zstd \
  gzip \
  -y

conda activate openspliceai_mouse
```

切换到你的分支并安装当前代码：

```bash
git clone https://github.com/Ethan-Fu-2000/OpenSpliceAI.git
cd OpenSpliceAI

git switch feat/imgt-splice-annotations

pip install -e .
```

确认：

```bash
openspliceai --help
python -c "import torch; print(torch.__version__); print(torch.cuda.is_available())"
```

公开 Release 通常可以直接下载。如果 `gh` 要求登录：

```bash
gh auth login
```

---

## 4. 只下载并解压数据

### 自动下载和解压

在仓库根目录执行：

```bash
bash scripts/download_mouse_immune_release.sh data
```

脚本自动完成：

1. 查找最新的 `mouse-immune-data-*` Release；
2. 下载所有 `.part` 分卷；
3. 校验每个分卷；
4. 合并为完整 `.tar.zst`；
5. 再校验完整压缩包；
6. 解压到 `data/mouse_immune/`。

解压后的关键结构：

```text
data/mouse_immune/
├── manifest.json
├── ncbi/
│   └── GCF_000001635.27/
│       └── ncbi_dataset/
│           └── data/
│               └── GCF_000001635.27/
│                   ├── *_genomic.fna
│                   ├── *_genomic.gff
│                   └── *_genomic.gbff
└── imgt/
    ├── gene_db/
    └── ligm_db/
        └── imgt.dat.Z
```

### 已经手工下载分卷时怎样解压

假设所有附件都在当前目录：

```bash
sha256sum -c PARTS_SHA256SUMS

cat mouse-immune-data.tar.zst.*.part \
  > mouse-immune-data.tar.zst

sha256sum -c ARCHIVE_SHA256SUM

mkdir -p data

tar --use-compress-program=zstd \
  -xf mouse-immune-data.tar.zst \
  -C data
```

确认：

```bash
test -f data/mouse_immune/manifest.json && echo "数据完整"
```

---

## 5. 最推荐：一键准备数据，先不训练

第一次建议只生成数据并检查，不要直接开始长时间训练：

```bash
python scripts/run_mouse_imgt_retraining.py \
  --skip-train \
  --flanking-size 2000 \
  --random-seed 42
```

这个命令会依次完成：

1. 本地没有数据时，从 Release 下载并解压；
2. 自动读取 `manifest.json`，定位 FASTA、GFF3 和 IMGT 文件；
3. 生成作者原始 NCBI protein-coding 基线；
4. 从 NCBI GFF3 提取 IG/TR 免疫位点；
5. 判断这些 NCBI 位点是否已经真正存在于基线 HDF5；
6. 将缺失的 NCBI 位点补入 train / validation；
7. 读取 IMGT functional 明确剪接位点；
8. 将 NCBI 未覆盖的 IMGT 位点补入 train / validation；
9. 原样复制作者 NCBI test；
10. 生成最终 `dataset_*.h5`；
11. 检查记录数、donor/acceptor 数量和 HDF5 shard；
12. 写出中文可读的审计与汇总 JSON。

主要输出：

```text
work/mouse_imgt_retraining/
├── run_config.json
├── ncbi_baseline/
│   ├── datafile_train.h5
│   ├── datafile_validation.h5
│   ├── datafile_test.h5
│   ├── dataset_train.h5
│   ├── dataset_validation.h5
│   └── dataset_test.h5
└── ncbi_imgt/
    ├── ncbi_imgt_splice_audit.json
    ├── training_data_summary.json
    ├── ncbi_supplement_train.h5
    ├── ncbi_supplement_validation.h5
    ├── imgt_supplement_train.h5
    ├── imgt_supplement_validation.h5
    ├── immune_supplement_train.h5
    ├── immune_supplement_validation.h5
    ├── datafile_train.h5
    ├── datafile_validation.h5
    ├── datafile_test.h5
    ├── dataset_train.h5
    ├── dataset_validation.h5
    └── dataset_test.h5
```

---

## 6. 训练前必须检查什么

查看总览：

```bash
cat work/mouse_imgt_retraining/ncbi_imgt/training_data_summary.json
```

查看详细审计：

```bash
less work/mouse_imgt_retraining/ncbi_imgt/ncbi_imgt_splice_audit.json
```

重点字段：

```text
supplemented_from_ncbi
imgt.supplemented_from_imgt
ncbi_already_in_baseline
policy.test_set
```

应该满足：

- `supplemented_from_ncbi` 可以为 0 或正数，取决于当前 NCBI 注释和作者基线；
- `imgt.supplemented_from_imgt` 应反映 NCBI 未覆盖的 IMGT 明确位点；
- `policy.test_set` 明确说明 test 没有加入任何免疫补充记录；
- donor motif 通常为 `GT`，acceptor motif 通常为 `AG`；
- 非典型 motif 只有在 NCBI 或 IMGT 明确标注时才保留；
- `IGHJ` donor 和 `IGHC` 首个 acceptor 应在审计明细中出现。

快速搜索：

```bash
grep -i -E '"gene": "IGHJ|"gene": "IGHM|"gene": "IGHG|"gene": "IGHA' \
  work/mouse_imgt_retraining/ncbi_imgt/ncbi_imgt_splice_audit.json \
  | head -50
```

---

## 7. 开始重新训练

### 8 GB 显存的稳妥起点

你的 RTX 5060 Ti 是 8 GB，建议先从：

- `flanking-size=2000`
- `batch-size=2`
- `num-gpus=1`
- `epochs=20`

开始：

```bash
python scripts/run_mouse_imgt_retraining.py \
  --flanking-size 2000 \
  --batch-size 2 \
  --num-gpus 1 \
  --epochs 20 \
  --random-seed 42 \
  --project-name mouse_ncbi_imgt \
  --exp-num 1
```

脚本发现数据已经生成后，不会重复处理，会直接检查并训练。

### 后台运行

```bash
mkdir -p logs

nohup python scripts/run_mouse_imgt_retraining.py \
  --flanking-size 2000 \
  --batch-size 2 \
  --num-gpus 1 \
  --epochs 20 \
  --random-seed 42 \
  --project-name mouse_ncbi_imgt \
  --exp-num 1 \
  > logs/mouse_ncbi_imgt_train.log 2>&1 &

echo $! > logs/mouse_ncbi_imgt_train.pid
```

查看日志：

```bash
tail -f logs/mouse_ncbi_imgt_train.log
```

查看 GPU：

```bash
watch -n 2 nvidia-smi
```

### 显存不足

出现 `CUDA out of memory` 时：

```bash
python scripts/run_mouse_imgt_retraining.py \
  --flanking-size 2000 \
  --batch-size 1 \
  --num-gpus 1 \
  --epochs 20
```

不要一开始使用 `flanking-size=10000`。它的显存、运行时间和数据读取开销都会显著增加。

---

## 8. 模型输出在哪里

默认输出：

```text
models/mouse_ncbi_imgt/
└── SpliceAI_mouse_ncbi_imgt_2000_1_rs42/
    └── 1/
        ├── models/
        │   ├── model_*.pt
        │   └── model_best.pt
        └── LOG/
            ├── TRAIN/
            ├── VAL/
            └── TEST/
```

最重要的模型通常是：

```text
models/mouse_ncbi_imgt/SpliceAI_mouse_ncbi_imgt_2000_1_rs42/1/models/model_best.pt
```

检查：

```bash
find models/mouse_ncbi_imgt -name 'model_best.pt' -o -name 'model_*.pt'
```

---

## 9. 使用新模型预测

假设模型为：

```bash
MODEL="models/mouse_ncbi_imgt/SpliceAI_mouse_ncbi_imgt_2000_1_rs42/1/models/model_best.pt"
```

对你的 5000 bp 小鼠序列预测：

```bash
openspliceai predict \
  --input-sequence your_mouse_sequence.fa \
  --model "${MODEL}" \
  --flanking-size 2000 \
  --output-dir predict_mouse_ncbi_imgt \
  --threshold 0.01
```

训练和预测的 `--flanking-size` 必须一致。训练使用 2000，预测也必须使用 2000。

---

## 10. 分步骤执行方式

不使用一键脚本时，可以手动执行。

### 10.1 下载解压

```bash
bash scripts/download_mouse_immune_release.sh data
```

### 10.2 从 manifest 获取路径

```bash
python - <<'PY'
import json
from pathlib import Path

root = Path("data/mouse_immune")
manifest = json.loads((root / "manifest.json").read_text())

print("FASTA:", root / manifest["ncbi"]["genome_fasta"])
print("GFF3 :", root / manifest["ncbi"]["annotation_gff3"])
print("IMGT :", root / manifest["imgt"]["flat_file"])
PY
```

### 10.3 生成作者 NCBI 基线

推荐仍使用一键脚本负责这一阶段，因为它会固定随机种子：

```bash
python scripts/run_mouse_imgt_retraining.py \
  --skip-train
```

### 10.4 单独执行 NCBI/IMGT 补充

将下面三个路径替换为 manifest 中的实际路径：

```bash
python scripts/prepare_mouse_immune_training.py \
  --ncbi-gff <GRCm39_GFF3> \
  --ncbi-fasta <GRCm39_FASTA> \
  --imgt-flat data/mouse_immune/imgt/ligm_db/imgt.dat.Z \
  --base-data-dir work/mouse_imgt_retraining/ncbi_baseline \
  --output-dir work/mouse_imgt_retraining/ncbi_imgt \
  --training-context-radius 5000 \
  --match-radius 20 \
  --identity-threshold 0.90 \
  --validation-ratio 0.10
```

---

## 11. 常用参数

| 参数 | 默认值 | 含义 |
|---|---:|---|
| `--flanking-size` | 2000 | 模型上下文长度，训练与预测必须一致 |
| `--batch-size` | 2 | 8 GB 显存稳妥起点 |
| `--num-gpus` | 1 | 当前训练使用的 GPU 数 |
| `--epochs` | 20 | 训练轮数 |
| `--random-seed` | 42 | 数据拆分和训练随机种子 |
| `--skip-train` | 关闭 | 开启后只准备数据 |
| `--force-baseline` | 关闭 | 删除并重建 NCBI 基线 |
| `--force-augment` | 关闭 | 删除并重建增强数据 |
| `--include-nonfunctional` | 关闭 | 是否加入 ORF/假基因；通常不要开启 |
| `--release-tag` | 最新 | 指定某一个云端数据 Release |

固定使用当前 Release：

```bash
python scripts/run_mouse_imgt_retraining.py \
  --release-tag mouse-immune-data-20260806-31072225004 \
  --skip-train
```

---

## 12. 重要限制

1. **IMGT 补充不是“凭距离猜 J-C 配对”**  
   只接受 NCBI 或 IMGT 明确标注的 splice site。

2. **模型学的是位点，不是重排组合类别**  
   它会学习 IGHJ donor、IGHC acceptor 周围的序列特征，但不会直接输出“该 J 与哪个 C 配对”。

3. **IMGT 序列与 GRCm39 可能不是完全同一等位基因背景**  
   代码使用 gene、位点类型和中心上下文比对来避免重复，并在审计 JSON 中保留来源。

4. **test 没有加入免疫补充数据**  
   这样能防止训练样本直接进入测试集。但若要专门评估 IG/TR 位点，还应该另建一个独立、按基因或等位基因隔离的免疫测试集。

5. **完整重训练可能耗时很长**  
   首次必须先用 `--skip-train` 检查数据和审计报告。

---

## 13. 最简操作顺序

```bash
# 1. 拉取代码
git fetch origin
git switch feat/imgt-splice-annotations
git pull

# 2. 安装
conda activate openspliceai_mouse
pip install -e .

# 3. 下载、解压、生成并检查训练数据
python scripts/run_mouse_imgt_retraining.py \
  --skip-train \
  --flanking-size 2000

# 4. 检查审计报告
cat work/mouse_imgt_retraining/ncbi_imgt/training_data_summary.json

# 5. 正式训练
mkdir -p logs
nohup python scripts/run_mouse_imgt_retraining.py \
  --flanking-size 2000 \
  --batch-size 2 \
  --num-gpus 1 \
  --epochs 20 \
  > logs/mouse_ncbi_imgt_train.log 2>&1 &
```
