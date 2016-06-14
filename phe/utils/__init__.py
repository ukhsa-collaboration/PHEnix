from math import floor

import psutil


def calculate_memory_for_sort():
    """Calculate available memory for ``samtools sort`` function.
    If there is enough memory, no temp files are created. **Enough**
    is defined as at least 1G per CPU.

    Returns
    -------
    sort_memory: str or None
        String to use directly with *-m* option in sort, or None.
    """
    avail_memory = psutil.virtual_memory().total
    avail_cpu = psutil.cpu_count()

    sort_memory = avail_memory / avail_cpu / 1024 ** 2

    # samtools default from documentation is 768M, to be conservative
    #    only use -m option when there is more than 1G per CPU.
    if sort_memory < 1024:
        sort_memory = None

    else:
        sort_memory = "%sG" % (int(floor(sort_memory / 1024)))

    return sort_memory


if __name__ == "__main__":
    print calculate_memory_for_sort()
