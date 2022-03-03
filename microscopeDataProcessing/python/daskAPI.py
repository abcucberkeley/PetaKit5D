# from dask_jobqueue import SLURMCluster
# from dask.distributed import Client, LocalCluster
import dask.array as da
import time

def daskZarrMaxProjection(zarrFullname, axis=2, log_directory=None):
    # cluster = SLURMCluster(queue="abc", project="co_abc", cores=24,memory="500GB",job_extra=["--qos=abc_normal"], dashboard_address=":8797", log_directory=log_directory)
    # client = Client(abccluster)
    # abccluster.scale(1)
    # cluster = LocalCluster(local_directory='/tmp/dask')
    # client = Client(cluster, timeout="50s")

    img = da.from_zarr(zarrFullname)
    t = time.time()
    MIP_z = img.max(axis=int(axis)).compute()
    elapsed = time.time() - t
    print("Elasped time: {} s".format(elapsed))
    # cluster.close()
    return MIP_z


def daskZarrPadArray(zarrFullname, OutputFullname, pad_width, mode='constant', log_directory=None):
    # cluster = SLURMCluster(queue="abc", project="co_abc", cores=24,memory="500GB",job_extra=["--qos=abc_normal"], dashboard_address=":8797", log_directory=log_directory)
    # client = Client(abccluster)
    # abccluster.scale(1)
    # cluster = LocalCluster(local_dir='/tmp/dask')
    # client = Client(cluster, timeout="50s")

    img = da.from_zarr(zarrFullname)
    t = time.time()
    da.pad(img, pad_width, mode='constant', constant_values=0).rechunk(img.chunksize).to_zarr(OutputFullname)
    elapsed = time.time() - t
    print("Elasped time: {} s".format(elapsed))
    # cluster.close()
    return 0

def main(zarrFullname=None, OutputFullname=None, pad_width=None, mode=None):

    if len(sys.argv) > 1:
        zarrFullname = sys.argv[1]
        OutputFullname = sys.argv[2]
        pad_width = sys.argv[3]

    print(sys.argv)
    print(zarrFullname)
    print(OutputFullname)
    print(pad_width)
    daskZarrPadArray(zarrFullname, outputFullpath, pad_width)

if __name__ == "__main__":
    main()

