import os
import random
from itertools import count
import hail as hl
from uuid import uuid4


wd = '/net/archive/groups/plggneuromol/ifpan-gosborcz-ukb/'

spark_master_host = os.environ.get('SPARK_MASTER_HOST')
spark_master_port = os.environ.get('SPARK_MASTER_PORT')
localfs_path = os.environ.get('SCRATCH_LOCAL') + '/'
scratch_path = os.environ.get('SCRATCH') + '/'
hail_log_uuid = str(uuid4())


def hl_init(
        master=f'spark://{spark_master_host}:{spark_master_port}',
        tmp_dir=os.path.join(scratch_path, 'hail-tmpdir'),
        default_reference='GRCh38',
        spark_conf={
            'spark.driver.memory': '40G',
            'spark.executor.memory': '80G',
            'spark.rpc.message.maxSize': '256',
        },
        log=os.path.join(scratch_path, f'slurm-log/hail-{hail_log_uuid}.log'),
        local_tmpdir=os.path.join(localfs_path, 'hail-local-tmpdir')):
    hl.init(
        #master=master,
        tmp_dir=tmp_dir,
        default_reference=default_reference,
        spark_conf=spark_conf,
        log=log,
        local_tmpdir=local_tmpdir,
    )


def tmpdir_path_iter(prefix=None):
    if int(os.getenv('SLURM_NNODES')) > 1 or os.getenv('HAIL_CHECKPOINT_ENV'):
        tmp_path = os.path.join(scratch_path, 'tmp/')
    else:
        tmp_path = os.path.join(localfs_path, 'tmp/')
    if prefix is None:
        prefix = f"{random.randrange(16 ** 4):04x}"
    counter = count()
    while True:
        path = os.path.join(
            tmp_path,
            f"{prefix}-{os.getenv('SLURM_JOBID')}-{str(next(counter))}"
        )
        if os.getenv('HAIL_CHECKPOINT_ENV'):
            print(f'HAIL_CHECKPOINT_PATH: {path}', flush=True)
        yield path
