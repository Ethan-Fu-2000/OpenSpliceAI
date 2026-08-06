"""
Filename: train.py
Author: Kuan-Hao Chao
Date: 2025-03-20
Description: Train the OpenSpliceAI model.
"""

import sys
import numpy as np
import torch
import torch.optim as optim
from openspliceai.train_base.openspliceai import *
from openspliceai.train_base.utils import *
from openspliceai.constants import *


def initialize_model_and_optim(
    device,
    flanking_size,
    epochs,
    scheduler,
    batch_size=None,
    num_gpus=None,
):
    """初始化模型，并允许外部安全地覆盖 batch size / GPU 数量。

    原实现把 ``N_GPUS`` 固定为 2，即使机器只有一张显卡也会把默认 batch
    size 乘以 2。这里改为默认读取实际 CUDA GPU 数；CPU/MPS 按 1 处理。
    ``batch_size`` 显式传入时优先使用，便于 8 GB 显存从较小值开始。
    """
    L = 32
    detected_gpus = torch.cuda.device_count() if torch.cuda.is_available() else 1
    N_GPUS = max(1, int(num_gpus or detected_gpus))

    W = np.asarray([11, 11, 11, 11])
    AR = np.asarray([1, 1, 1, 1])
    base_batch_size = 18
    if int(flanking_size) == 80:
        W = np.asarray([11, 11, 11, 11])
        AR = np.asarray([1, 1, 1, 1])
        base_batch_size = 18
    elif int(flanking_size) == 400:
        W = np.asarray([11, 11, 11, 11, 11, 11, 11, 11])
        AR = np.asarray([1, 1, 1, 1, 4, 4, 4, 4])
        base_batch_size = 18
    elif int(flanking_size) == 2000:
        W = np.asarray([
            11, 11, 11, 11, 11, 11, 11, 11,
            21, 21, 21, 21,
        ])
        AR = np.asarray([
            1, 1, 1, 1, 4, 4, 4, 4,
            10, 10, 10, 10,
        ])
        base_batch_size = 12
    elif int(flanking_size) == 10000:
        W = np.asarray([
            11, 11, 11, 11, 11, 11, 11, 11,
            21, 21, 21, 21, 41, 41, 41, 41,
        ])
        AR = np.asarray([
            1, 1, 1, 1, 4, 4, 4, 4,
            10, 10, 10, 10, 25, 25, 25, 25,
        ])
        base_batch_size = 6
    else:
        raise ValueError("flanking_size 必须是 80、400、2000 或 10000")

    BATCH_SIZE = int(batch_size) if batch_size else base_batch_size * N_GPUS
    if BATCH_SIZE < 1:
        raise ValueError("batch_size 必须大于 0")

    CL = 2 * np.sum(AR * (W - 1))
    print("\033[1mContext nucleotides: %d\033[0m" % CL)
    print("\033[1mSequence length (output): %d\033[0m" % SL)
    print(
        f"\033[1mTraining batch size: {BATCH_SIZE}; "
        f"GPU count used for clipping: {N_GPUS}\033[0m"
    )

    model = SpliceAI(L, W, AR).to(device)
    print(model, file=sys.stderr)
    optimizer = torch.optim.AdamW(model.parameters(), lr=1e-3)
    if scheduler == "MultiStepLR":
        milestones = sorted(
            {
                max(1, epochs - 5),
                max(1, epochs - 4),
                max(1, epochs - 3),
                max(1, epochs - 2),
                max(1, epochs - 1),
            }
        )
        scheduler_obj = optim.lr_scheduler.MultiStepLR(
            optimizer,
            milestones=milestones,
            gamma=0.5,
        )
    elif scheduler == "CosineAnnealingWarmRestarts":
        scheduler_obj = torch.optim.lr_scheduler.CosineAnnealingWarmRestarts(
            optimizer,
            T_0=5,
            T_mult=1,
            eta_min=1e-5,
            last_epoch=-1,
        )
    else:
        raise ValueError(f"不支持的 scheduler：{scheduler}")

    params = {
        "L": L,
        "W": W,
        "AR": AR,
        "CL": CL,
        "SL": SL,
        "BATCH_SIZE": BATCH_SIZE,
        "N_GPUS": N_GPUS,
    }
    return model, optimizer, scheduler_obj, params


def train(args):
    """从头训练 SpliceAI 模型。"""
    print("Running OpenSpliceAI with 'train' mode")
    device = setup_environment(args)
    (
        model_output_base,
        log_output_train_base,
        log_output_val_base,
        log_output_test_base,
    ) = initialize_paths(args)
    train_h5f, valid_h5f, test_h5f, batch_num = load_datasets(args)
    train_idxs, val_idxs, test_idxs = generate_indices(
        train_h5f,
        valid_h5f,
        test_h5f,
    )

    model, optimizer, scheduler, params = initialize_model_and_optim(
        device,
        args.flanking_size,
        args.epochs,
        args.scheduler,
        batch_size=getattr(args, "batch_size", None),
        num_gpus=getattr(args, "num_gpus", None),
    )
    params["RANDOM_SEED"] = args.random_seed
    train_metric_files = create_metric_files(log_output_train_base)
    valid_metric_files = create_metric_files(log_output_val_base)
    test_metric_files = create_metric_files(log_output_test_base)
    train_model(
        model,
        optimizer,
        scheduler,
        train_h5f,
        valid_h5f,
        test_h5f,
        train_idxs,
        val_idxs,
        test_idxs,
        model_output_base,
        args,
        device,
        params,
        train_metric_files,
        valid_metric_files,
        test_metric_files,
    )
    train_h5f.close()
    valid_h5f.close()
    test_h5f.close()
