import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
import argparse

def main():
    """Convert native data for a single sampler point to a time series in ASCII."""
    parser = argparse.ArgumentParser(
        description="Convert native data for a single sampler point to a time series in ASCII"
    )
    parser.add_argument(
        "-ks",
        "--kynema_sgf_path",
        help="path to amr-wind directory",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-d",
        "--run_dir",
        help="Directory where run took place, where post-processing directory is",
        required=True,
        type=str,
    )
    parser.add_argument(
        "-p",
        "--post_proc_dir",
        help="Name of post-processing directory",
        default="post_processing",
        type=str,
    )
    parser.add_argument(
        "-l",
        "--label",
        help="Top-level sampling label name",
        default="sampling",
        type=str,
    )
    parser.add_argument(
        "-s",
        "--sampler",
        help="Bottom-level sampling label name",
        default="* use first *",
        type=str,
    )
    parser.add_argument(
        "-i",
        "--index",
        help="Particle index within specified sampler",
        default=0,
        type=int,
    )
    args = parser.parse_args()

    ########## Input arguments from command line ##############
    print("\nParameters to be used:")
    print("   run directory path: ", os.path.realpath(args.run_dir))
    print("   Kynema-SGF path: ", os.path.realpath(args.kynema_sgf_path))
    print("   post-processing directory name: ", args.post_proc_dir)
    print("   top-level sampling label name: ", args.label)
    ###########################################################

    sys.path.append(args.kynema_sgf_path+'/tools/')
    import amrex_particle
    from amrex_particle import AmrexParticleFile

    post_proc_dir = args.run_dir + "/" + args.post_proc_dir
    pfile = AmrexParticleFile(args.run_dir)

    os.chdir(post_proc_dir)
    contents = glob.glob(args.label + "*")
    directories = [path for path in contents if os.path.isdir(path)]
    sdirs = sorted(directories)
    n_dir = len(sdirs)
    if (n_dir == 0):
        print("ERROR: No matching sampling directories found, exiting!")
        sys.exit(1)

    # Separate by length, then order again (because 10000 would be next to 100000)
    max_len = len(sdirs[0])
    min_len = max_len
    for n in range(len(sdirs)):
        max_len = max(max_len,len(sdirs[n]))
        min_len = min(min_len,len(sdirs[n]))

    sdirs_init = list(sdirs)
    if (max_len != min_len):
        cur_len = min_len
        idx = 0
        for inc_len in range(max_len - min_len + 1):
            cur_len += inc_len
            for n in range(len(sdirs_init)):
                if (len(sdirs_init[n]) == cur_len):
                    sdirs[idx] = sdirs_init[n]
                    idx += 1

    max_step = int(sdirs[n_dir-1][len(args.label):len(sdirs[n_dir-1])])

    # Get stuff for header and array sizes using first available directory
    n0 = 0
    step = int(sdirs[n0][len(args.label):len(sdirs[n0])])
    pt = pfile.load(step, label=args.label, root_dir = ".")
    pt.parse_header()
    pt.parse_info()
    output = np.zeros((n_dir, 1 + 3 + pt.num_reals + pt.num_ints - 3))
    v_str = "x y z"
    for rv in range(pt.num_reals):
        v_str += " " + pt.real_var_names[rv]
    for iv in range(pt.num_ints - 3):
        v_str += " " + pt.int_var_names[iv + 3]
    header = "t " + v_str
    samplers = pt.info["samplers"]
    s_index = 0
    if (args.sampler[0] != "*"):
        for s in range(len(samplers)):
            if (args.sampler == samplers[s]["label"]):
                s_index = s
    else:
        args.sampler = samplers[0]["label"]

    ########## Input arguments from command line ##############
    print("   bottom-level sampling label name: ", args.sampler)
    print("   particle index within sampler: ", args.index)
    ###########################################################

    fname = args.label + "_" + args.sampler + "_" + str(args.index) + "_time_series.txt"
    print("Output file name: ", fname)
    print("Output file name with path: ", os.path.realpath(fname))
    print()

    # Loop through files to get time series data
    for n in range(n_dir):
        step = int(sdirs[n][len(args.label):len(sdirs[n])])
        print(str(step) + " out of " + str(max_step))
        pt = pfile.load(step, label=args.label, root_dir = ".")
        pt.parse_header()
        pt.parse_info()
        pt.load_binary_data()
        data = pt.df
        
        output[n, 0] = pt.info["time"]
        
        for p in range(pt.num_particles):
            if (data.set_id[p] == s_index and data.probe_id[p] == args.index):
                output[n, 1] = data.xco[p]
                output[n, 2] = data.yco[p]
                output[n, 3] = data.zco[p]
                for rv in range(pt.num_reals):
                    output[n, 4+rv] = data[pt.real_var_names[rv]][p]
                for iv in range(pt.num_ints - 3):
                    output[n, 4+pt.num_reals+iv] = data[pt.int_var_names[iv+3]][p]
                    
    np.savetxt(fname, output[0:len(sdirs)],header=header)

if __name__ == "__main__":
    main()  