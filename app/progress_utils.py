# app/progress_utils.py

from tqdm import tqdm

def update_progress_bar(pbar, total_steps):
    progress = 100 / total_steps
    if pbar.n + progress > pbar.total:
        pbar.update(pbar.total - pbar.n)  # Final clamp to ensure no overflow
    else:
        pbar.update(progress)


"""
Copyright (c) 2024 VECTOR FUTURES LTD
All rights reserved.
This file is part of the Vector Futures Codon Optimization App and may not be copied, distributed, or modified without express written permission.
"""