import numpy as np
import matplotlib.pyplot as plt
import sys
import glob
import os
import argparse

def main():
    """Convert native sampling data to a structured ASCII file for every sampler and output step"""
    parser = argparse.ArgumentParser(
        description="Convert native sampling data to a structured ASCII file for every sampler and output step"
    )
    parser.add_argument(
        "-ks",
        "--kynema_sgf_path",
        help="path to kynema-sgf directory",
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
        "-n",
        "--number_of_digits",
        help="Fixed number of digits denoting output step",
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

    # if fixed length is specified, prune directories that do not fit
    skip_flag = np.zeros(len(sdirs))
    ndir = len(sdirs)
    max_dir = ndir - 1
    if (args.number_of_digits > 0):
        ndir = 0
        max_dir = ndir - 1
        for n in range(len(sdirs)):
            if (int(len(sdirs[n])) != int(len(args.label) + args.number_of_digits)):
                skip_flag[n] = 1
            else:
                ndir += 1
                max_dir = n

    print("Number of folders to read: " + str(ndir))
    if (max_dir < 0):
        print("ERROR: No matching sampling directories found, exiting!")
        sys.exit(1)
    max_step = int(sdirs[max_dir][len(args.label):len(sdirs[max_dir])])

    samplers_printed = False
    for n in range(len(sdirs)):
        if (skip_flag[n] == 1):
            continue

        step = int(sdirs[n][len(args.label):len(sdirs[n])])
        if (samplers_printed):
            print(str(step) + " out of " + str(max_step))
        pt = pfile.load(step, label=args.label, root_dir = ".")
        pt.parse_header()
        pt.parse_info()
        pt.load_binary_data()
        data = pt.df
        samplers = pt.info["samplers"]
        if (not samplers_printed):
            samplers_printed = True
            print("Sampler labels found in first step:")
            for sampler in samplers:
                print("   " + sampler["label"])
            print("Reading folder " + str(step) + " out of " + str(max_step))
        
        t_str = str(pt.info["time"])
        v_str = "x y z"
        
        for rv in range(pt.num_reals):
            v_str += " " + pt.real_var_names[rv]
        for iv in range(pt.num_ints - 3):
            v_str += " " + pt.int_var_names[iv + 3]

        header = "t = " + t_str + "\n" + v_str
            
        # Have to go with max possible size, subsampler size not stored
        output = np.zeros((pt.num_particles, 3 + pt.num_reals + pt.num_ints - 3))
        
        for s in range(len(samplers)):
            np_s = 0
            for p in range(pt.num_particles):
                if data.set_id[p] == samplers[s]["index"]:
                    np_s += 1
                    pid = data.probe_id[p]
                    output[pid, 0] = data.xco[p]
                    output[pid, 1] = data.yco[p]
                    output[pid, 2] = data.zco[p]
                    for rv in range(pt.num_reals):
                        output[pid, 3+rv] = data[pt.real_var_names[rv]][p]
                    for iv in range(pt.num_ints - 3):
                        output[pid, 3+pt.num_reals+iv] = data[pt.int_var_names[iv+3]][p]
                    

            fname = sdirs[n] + "_" + samplers[s]["label"]+ ".txt"
            np.savetxt(fname,output[0:np_s,:],header=header)

if __name__ == "__main__":
    main() 
