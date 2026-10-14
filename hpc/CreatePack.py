import os

def generate_matlab_command(
    N, K, D, G, M, P_target, seed,
    plotit, x_mult, y_mult, z_mult,
    calc_eig, save_path
):
    plotit_str   = "true" if plotit else "false"
    calc_eig_str = "true" if calc_eig else "false"

    # z_mult decides 2D vs 3D inside MATLAB (boolThreeD = (z_mult ~= 0)).
    # pack now takes 13 positional args (no name-value 'hertzian' pair);
    # hertzian is only selectable via the optional trailing options struct,
    # which we don't pass here, so it defaults to false (Hooke).
    #
    # 3D repeat-tile requires at least one of x/y/z mult to differ from 1,
    # and pack's default options (13-arg call) include saveFullState = true,
    # so pack writes the pre-cleanRats '_Full' tile source and then calls
    # packRepeatTile to produce the 3D superlattice in one MATLAB invocation.
    return (
        f"matlab -nodisplay -nosplash -r \"addpath('./src/'); try; "
        f"pack({N}, {K}, {D}, {G}, {M}, {P_target}, {seed}, "
        f"{plotit_str}, {x_mult}, {y_mult}, {z_mult}, "
        f"{calc_eig_str}, '{save_path}'); "
        f"catch e; disp(e.message); end; exit\""
    )

def main():
    # Physical parameters (unchanged).
    K = 100        # spring constant
    D = 1          # small-particle diameter
    G = 1.4        # large/small diameter ratio
    M = 1          # particle mass
    P_target = 0.001   # target pressure -- matches the previously made Width20 packings
    seed = 1       # RNG seed
    plotit = False
    calc_eig = False
    save_path = "./data/packings/3d/hooke/"

    # 3D base packings to generate and tile.
    # Width in the filename = round(N^(1/3)); a "WxWxW" base is N = W^3.
    # The tiling multiplier x*y*z sets the final particle count:
    #   20^3 = 8000   -> tile 6 x 5 x 3 = 90   -> 720,000
    #   40^3 = 64000  -> tile 2 x 2 x 3 = 12   -> 768,000
    # (both "around N = 700,000")
    base_packings = [
        (8000, 6, 5, 3),   # 20x20x20 box, 3D, tiled 6x5x3 -> 720,000
        (64000, 2, 2, 3),  # 40x40x40 box, 3D, tiled 2x2x3 -> 768,000
    ]

    output_file = os.path.join(os.path.dirname(__file__), "commandsPack.txt")

    with open(output_file, "w") as file:
        for N, x_mult, y_mult, z_mult in base_packings:
            command = generate_matlab_command(
                    N, K, D, G, M,
                    P_target, seed,
                    plotit, x_mult, y_mult, z_mult,
                    calc_eig, save_path
                    )
            file.write(command + "\n")

    print(f"Commands written to: {output_file}")
    for N, x_mult, y_mult, z_mult in base_packings:
        width = int(round(N ** (1.0 / 3.0)))
        print(f"  N={N} ({width}^3) tile {x_mult}x{y_mult}x{z_mult} "
              f"(x{x_mult*y_mult*z_mult}) -> N_tiled={N * x_mult * y_mult * z_mult}")

if __name__ == "__main__":
    main()
