# 小鼠 OpenSpliceAI：NCBI + IMGT 严格重训练中文指南

## 1. 最终只使用两个来源

本流程只接入：

1. **NCBI RefSeq GRCm39**：`GCF_000001635.27`，作为主体训练集；
2. **IMGT**：只读取明确标注的 donor/acceptor，默认只用 functional 小鼠记录，并要求位点中心序列在 GRCm39 中全长零错配、唯一定位。

**GENCODE 不下载、不解析、不加入训练，也不作为补缺来源。**

同时不接入：

- mm10/GRCm38；
- liftOver 坐标；
- VastDB；
- recount3；
- ENCODE novel junction；
- 纯算法预测位点。

## 2. 数据优先级

```text
NCBI GRCm39 原始基线
    ↓
基线实际遗漏的 NCBI IG/TR 明确位点
    ↓
精确匹配 GRCm39、且 NCBI 未覆盖的 IMGT functional 位点
```

原始测试集保持不变。凡是映射到测试染色体的补充位点，一律排除。

## 3. 更新代码

```bash
cd ~/Documents/yanfa/OpenSpliceAI

git fetch origin
git switch feat/imgt-splice-annotations
git pull

conda activate openspliceai
pip install -e .
```

确认环境：

```bash
python -c "import torch; print(torch.__version__); print(torch.cuda.is_available())"
```

严格 IMGT 定位需要 `mappy`。如果环境里没有：

```bash
pip install mappy
```

## 4. 下载并解压纯 NCBI+IMGT 数据

一键流程会自动下载最新 Release。也可以单独执行：

```bash
bash scripts/download_mouse_immune_release.sh data
```

下载脚本会：

1. 下载所有分卷；
2. 校验每个分卷 SHA-256；
3. 合并完整压缩包；
4. 校验完整压缩包 SHA-256；
5. 删除本地旧版 `data/mouse_immune/`；
6. 解压新数据。

第 5 步会清除旧版残留的 `gencode/` 目录。

解压后结构应为：

```text
data/mouse_immune/
├── manifest.json
├── ncbi/
└── imgt/
```

不应再有：

```text
data/mouse_immune/gencode/
```

检查清单：

```bash
python - <<'PY'
import json
from pathlib import Path

m = json.loads(Path("data/mouse_immune/manifest.json").read_text())
print("profile:", m["policy"]["profile"])
print("sources:", m["policy"]["sources"])
print("gencode:", m["policy"]["gencode_included"])
print("NCBI:", m["ncbi"]["assembly_accession"], m["ncbi"]["assembly_name"])
print("IMGT:", m["imgt"]["gene_db_release"], m["imgt"]["ligm_db_release"])
PY
```

必须看到：

```text
profile: ncbi_imgt_strict_v2
sources: ['NCBI RefSeq', 'IMGT']
gencode: False
NCBI: GCF_000001635.27 GRCm39
```

## 5. 第一次只生成训练数据

```bash
python scripts/run_mouse_strict_retraining.py \
  --skip-train \
  --flanking-size 2000 \
  --random-seed 42
```

流程会自动完成：

1. 检查数据清单；
2. 生成作者原始 NCBI protein-coding 基线；
3. 检查 NCBI IG/TR 位点是否已进入基线；
4. 补入基线遗漏的 NCBI IG/TR 位点；
5. 读取 IMGT functional 小鼠明确位点；
6. 将位点中心序列零错配定位到 GRCm39；
7. 默认仅保留唯一定位位点；
8. 排除测试染色体补充位点；
9. 排除已被 NCBI 覆盖的 IMGT 位点；
10. 生成最终 HDF5 和审计报告。

结果目录：

```text
work/mouse_ncbi_imgt_strict_retraining/
├── ncbi_baseline/
└── ncbi_imgt_strict/
    ├── ncbi_imgt_strict_splice_audit.json
    ├── training_data_summary.json
    ├── datafile_train.h5
    ├── datafile_validation.h5
    ├── datafile_test.h5
    ├── dataset_train.h5
    ├── dataset_validation.h5
    └── dataset_test.h5
```

## 6. 查看审计结果

数据摘要：

```bash
cat \
  work/mouse_ncbi_imgt_strict_retraining/ncbi_imgt_strict/training_data_summary.json
```

详细报告：

```bash
less \
  work/mouse_ncbi_imgt_strict_retraining/ncbi_imgt_strict/ncbi_imgt_strict_splice_audit.json
```

重点字段：

```text
policy.gencode
policy.test_set

ncbi.supplemented
ncbi.removed_test_chromosome

imgt.genome_match.accepted_exact_unique_site_count
imgt.genome_match.rejected_no_exact_grcm39_match
imgt.genome_match.rejected_multiple_exact_grcm39_matches
imgt.removed_test_chromosome
imgt.duplicate_filter.supplemented_from_imgt
```

应确认：

```text
policy.gencode = not used
```

## 7. 正式训练

RTX 5060 Ti 8 GB 建议从 `batch-size 2` 开始：

```bash
mkdir -p logs

nohup python scripts/run_mouse_strict_retraining.py \
  --flanking-size 2000 \
  --batch-size 2 \
  --num-gpus 1 \
  --epochs 20 \
  --random-seed 42 \
  --project-name mouse_ncbi_imgt_strict \
  --exp-num 1 \
  > logs/mouse_ncbi_imgt_strict.log 2>&1 &

echo $! > logs/mouse_ncbi_imgt_strict.pid
```

查看日志：

```bash
tail -f logs/mouse_ncbi_imgt_strict.log
```

查看显卡：

```bash
watch -n 2 nvidia-smi
```

显存不足时改为：

```text
--batch-size 1
```

## 8. 默认不要开启的选项

下面两个参数默认关闭：

```text
--include-nonfunctional
--allow-nonunique-imgt
```

原因：

- `--include-nonfunctional` 会加入 ORF/假基因记录；
- `--allow-nonunique-imgt` 会允许多重定位，降低坐标可信度。

当前目标是少而精，因此不建议开启。

## 9. 模型输出

默认目录：

```text
models/mouse_ncbi_imgt_strict/
```

查找最佳模型：

```bash
find models/mouse_ncbi_imgt_strict \
  -name 'model_best.pt' -o -name 'model_*.pt'
```

## 10. 最终原则

本流程不是为了收集最多的可变剪切位点，而是为了得到一套可审计的免疫增强模型：

```text
NCBI 作为主体
IMGT 只补明确位点
必须严格对应 GRCm39
不使用 GENCODE
不修改测试集
```
