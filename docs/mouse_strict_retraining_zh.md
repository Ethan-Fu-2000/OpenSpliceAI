# 小鼠 OpenSpliceAI 严格重训练中文指南

## 1. 只使用三个来源

本流程只接入：

1. **NCBI RefSeq GRCm39**：`GCF_000001635.27`，作为主体训练集；
2. **GENCODE Mouse M39 Basic**：同为 `GRCm39`，只接受 Basic、level 1/2 和规范剪接 motif；
3. **IMGT**：只读取明确的 donor/acceptor 标注，默认只用 functional 小鼠记录，并要求位点中心序列在 GRCm39 中全长、零错配、唯一定位。

不接入旧版 mm10/GRCm38、liftOver 数据、VastDB、recount3、ENCODE novel junction 或纯算法预测位点。

## 2. 数据优先级

```text
NCBI GRCm39 基线
    ↓
基线实际缺失的 NCBI IG/TR 位点
    ↓
NCBI 未覆盖的 GENCODE M39 Basic 位点
    ↓
精确唯一匹配 GRCm39、且前两者未覆盖的 IMGT 位点
```

能定位到原始测试集染色体的补充位点全部排除，测试集保持原样。

## 3. 更新并安装代码

```bash
cd ~/Documents/yanfa/OpenSpliceAI

git fetch origin
git switch feat/imgt-splice-annotations
git pull

conda activate openspliceai
pip install -e .
```

确认：

```bash
python -c "import torch; print(torch.__version__); print(torch.cuda.is_available())"
```

## 4. 下载并解压云端数据

一键脚本会自动下载最新 Release。也可以先单独下载：

```bash
bash scripts/download_mouse_immune_release.sh data
```

解压后应看到：

```text
data/mouse_immune/
├── manifest.json
├── ncbi/
├── gencode/
│   └── gencode.vM39.basic.annotation.gtf.gz
└── imgt/
```

检查版本：

```bash
python - <<'PY'
import json
from pathlib import Path

m = json.loads(Path("data/mouse_immune/manifest.json").read_text())
print("NCBI:", m["ncbi"]["assembly_accession"], m["ncbi"]["assembly_name"])
print("GENCODE:", m["gencode"]["release"], m["gencode"]["assembly"])
print("GENCODE file:", m["gencode"]["basic_annotation_gtf"])
print("IMGT:", m["imgt"]["gene_db_release"], m["imgt"]["ligm_db_release"])
PY
```

必须看到：

```text
NCBI: GCF_000001635.27 GRCm39
GENCODE: M39 GRCm39
```

## 5. 第一次只准备数据

```bash
python scripts/run_mouse_strict_retraining.py \
  --skip-train \
  --flanking-size 2000 \
  --random-seed 42
```

主要结果：

```text
work/mouse_gencode_imgt_retraining/
├── ncbi_baseline/
└── ncbi_gencode_imgt/
    ├── ncbi_gencode_imgt_splice_audit.json
    ├── training_data_summary.json
    ├── datafile_train.h5
    ├── datafile_validation.h5
    ├── datafile_test.h5
    ├── dataset_train.h5
    ├── dataset_validation.h5
    └── dataset_test.h5
```

## 6. 检查审计结果

```bash
cat work/mouse_gencode_imgt_retraining/ncbi_gencode_imgt/training_data_summary.json
```

详细报告：

```bash
less work/mouse_gencode_imgt_retraining/ncbi_gencode_imgt/ncbi_gencode_imgt_splice_audit.json
```

重点字段：

```text
policy.assembly
policy.test_set

ncbi.supplemented
ncbi.removed_test_chromosome

gencode.parse.header_verified_release
gencode.parse.header_verified_assembly
gencode.parse.accepted_transcripts
gencode.parse.rejected_noncanonical_intron_pairs
gencode.supplemented
gencode.removed_test_chromosome

imgt.genome_match.accepted_exact_unique_site_count
imgt.genome_match.rejected_no_exact_grcm39_match
imgt.genome_match.rejected_multiple_exact_grcm39_matches
imgt.removed_test_chromosome
imgt.duplicate_filter.supplemented_from_imgt
```

## 7. 正式训练

RTX 5060 Ti 8 GB 建议先用：

```bash
mkdir -p logs

nohup python scripts/run_mouse_strict_retraining.py \
  --flanking-size 2000 \
  --batch-size 2 \
  --num-gpus 1 \
  --epochs 20 \
  --random-seed 42 \
  --project-name mouse_ncbi_gencode_imgt \
  --exp-num 1 \
  > logs/mouse_ncbi_gencode_imgt.log 2>&1 &

echo $! > logs/mouse_ncbi_gencode_imgt.pid
```

查看：

```bash
tail -f logs/mouse_ncbi_gencode_imgt.log
watch -n 2 nvidia-smi
```

显存不足时改为：

```bash
--batch-size 1
```

## 8. 强制重新生成

```bash
python scripts/run_mouse_strict_retraining.py \
  --skip-train \
  --force-baseline \
  --force-augment \
  --flanking-size 2000
```

只重做增强数据：

```bash
python scripts/run_mouse_strict_retraining.py \
  --skip-train \
  --force-augment \
  --flanking-size 2000
```

## 9. 使用新模型预测

查找模型：

```bash
find models/mouse_ncbi_gencode_imgt \
  -name 'model_best.pt' -o -name 'model_*.pt'
```

预测：

```bash
MODEL="models/mouse_ncbi_gencode_imgt/SpliceAI_mouse_ncbi_gencode_imgt_2000_1_rs42/1/models/model_best.pt"

openspliceai predict \
  --input-sequence your_mouse_sequence.fa \
  --model "${MODEL}" \
  --flanking-size 2000 \
  --output-dir predict_mouse_strict \
  --threshold 0.01
```

训练和预测的 `--flanking-size` 必须一致。
