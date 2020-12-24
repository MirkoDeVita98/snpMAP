import os


def prepare_directiories(cache_dir, results_dir):
    if not os.path.exists(cache_dir):
        os.mkdir(cache_dir)
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)


def clear_cache(mode, cache_dir):
    content = os.listdir(cache_dir)
    for fname in content:
        if mode in fname:
            os.remove(os.path.join(cache_dir, fname))