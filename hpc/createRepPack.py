import os
import argparse


def generate_matlab_command(
    N, K, P_target, scalWidthFactor, seed,
    scalXMult, scalYMult, calc_eig,
    strInPath, strSavePath,
    scalZMult, boolHertzian
):
    """Build one `matlab -r` command line that runs src/packRepeatTile.m.

    packRepeatTile loads a saved base packing and repeats (tiles) it in x, y
    and/or z to build a superlattice. The 2D vs 3D path is selected the same
    way pack.m selects it: scalZMult ~= 0 selects 3D, scalZMult == 0 selects
    2D. For 3D the base tile source must be the full (pre-cleanRats) file that
    pack.m writes when options.saveFullState = true (3D_N<N>_P<P>_Width<W>_
    Seed<s>_Full[_Hertz].mat); packRepeatTile3D loads that file and saves the
    superlattice under the normal packing name.

    calc_eig_str must be a MATLAB boolean ("true"/"false"); strInPath and
    strSavePath are the directory of the base tile / the save location.
    """
    plotit_str   = "false"  # packRepeatTile always opens a figure for the
                            # result; pass calc_eig_str as the 8th positional
                            # arg (the calc_eig flag), not a separate plot flag.
    calc_eig_str = "true" if calc_eig else "false"
    z_str        = str(scalZMult)
    hertz_str    = "true" if boolHertzian else "false"

    return (
        f"matlab -nodisplay -nosplash -r \"addpath('./src/'); try; "
        f"packRepeatTile({N}, {K}, {P_target}, {scalWidthFactor}, {seed}, "
        f"{scalXMult}, {scalYMult}, {calc_eig_str}, '{strInPath}', '{strSavePath}', "
        f"{z_str}, {hertz_str}); "
        f"catch e; disp(e.message); end; exit\""
    )


def main():
    # Parameter grid. scalZMult selects the path:
    #   scalZMult = 0  -> 2D packing (tiles x/y only)
    #   scalZMult ~= 0 -> 3D packing (tiles x/y/z; needs the _Full source file)
    N_values         = [1000]
    K_values         = [100]
    P_target_values = [.0001, .0003, .001, .003, .01, .03, .1]
    scalWidth_values = [10]     # 3D: round(N^(1/3)); 2D: round(sqrt(N))
    seed_values      = [1,2]

    # Repeat multipliers (how many times each tile is copied along each axis).
    scalXMult_values = [640]
    scalYMult_values = [1]
    scalZMult_values = [1]      # 0 -> 2D, >=1 -> 3D

    calc_eig     = False
    boolHertzian = False

    strInPath    = "./data/packings/3d/hooke/"
    strSavePath  = "./data/packings/3d/hooke/"

    output_file = os.path.join(os.path.dirname(__file__), "commandsRepPack.txt")

    with open(output_file, "w") as file:
        for N in N_values:
            for K in K_values:
                for P_target in P_target_values:
                    for scalWidth in scalWidth_values:
                        for seed in seed_values:
                            for scalXMult in scalXMult_values:
                                for scalYMult in scalYMult_values:
                                    for scalZMult in scalZMult_values:
                                        command = generate_matlab_command(
                                            N, K, P_target, scalWidth, seed,
                                            scalXMult, scalYMult, calc_eig,
                                            strInPath, strSavePath,
                                            scalZMult, boolHertzian
                                        )
                                        file.write(command + "\n")

    print(f"Commands written to: {output_file}")


if __name__ == "__main__":
    main()
