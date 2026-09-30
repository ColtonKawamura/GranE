import os

def generate_matlab_command(
    K, M, Bv, w_D, N, P, W, seed,
    in_path, out_path,
    shear, fullSpectrum, maxAmpTracking, cleanRats
):
    """Build one `matlab -r` command line that runs src/simMD.m for a single
    (packing, parameter) combination.

    simMD's options argument is a positional struct (see src/simMD.m /
    src/pack.m), so struct(...) is passed as the 11th positional argument.
    Only the requested flags are set; simMD's arguments block backfills the
    rest to defaults. An empty struct() is valid when no flags are set.
    """
    fields = []
    if shear:
        fields.append("'shear', true")
    if fullSpectrum:
        fields.append("'fullSpectrum', true")
    if maxAmpTracking:
        fields.append("'maxAmpTracking', true")
    if cleanRats:
        fields.append("'cleanRats', true")
    options_str = "struct(" + ", ".join(fields) + ")"

    return (
        f"matlab -nodisplay -nosplash -r \"addpath('./src/'); try; "
        f"simMD({K}, {M}, {Bv}, {w_D}, {N}, {P}, {W}, {seed}, "
        f"'{in_path}', '{out_path}', {options_str}); "
        f"catch e; disp(getReport(e)); end; exit\""
    )

def main():
    # Parameter grids (mirror CreatePack.py: one .mat file per combination)
    K_values        = [100]
    M_values        = [1]
    Bv_values       = [1, 0.3, 0.1, 0.03, 0.01]
    w_D_values      = [10, 3, 1, .3, .1]
    N_values        = [160000]
    P_values        = [0.1, 0.03, 0.01, 0.003, 0.001, 0.0003, 0.0001]
    W_values        = [10]
    seed_values     = [1]

    in_path  = "./data/packings/3d/hooke/"
    out_path = "./data/simMD/3d/hooke/"

    shear          = False
    fullSpectrum   = False
    maxAmpTracking = False
    cleanRats      = False

    output_file = os.path.join(os.path.dirname(__file__), "commandsSim.txt")

    with open(output_file, "w") as file:
        for K in K_values:
            for M in M_values:
                for Bv in Bv_values:
                    for w_D in w_D_values:
                        for N in N_values:
                            for P in P_values:
                                for W in W_values:
                                    for seed in seed_values:
                                        command = generate_matlab_command(
                                            K, M, Bv, w_D, N, P, W, seed,
                                            in_path, out_path,
                                            shear, fullSpectrum, maxAmpTracking, cleanRats
                                        )
                                        file.write(command + "\n")

    print(f"Commands written to: {output_file}")

if __name__ == "__main__":
    main()
