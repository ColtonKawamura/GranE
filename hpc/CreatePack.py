import numpy as np

def generate_matlab_command(
    N, K, D, G, M, P_target, seed,
    plotit, x_mult, y_mult, z_mult,
    calc_eig, save_path, hertzian
):
    plotit_str   = "true" if plotit else "false"
    calc_eig_str = "true" if calc_eig else "false"
    hertzian_str = "true" if hertzian else "false"

    # z_mult decides 2D vs 3D inside MATLAB (boolThreeD = (z_mult ~= 0))
    return (
        f"matlab -nodisplay -nosplash -r \"addpath('./src/'); try; "
        f"pack({N}, {K}, {D}, {G}, {M}, {P_target}, {seed}, "
        f"{plotit_str}, {x_mult}, {y_mult}, {z_mult}, "
        f"{calc_eig_str}, '{save_path}', 'hertzian', {hertzian_str}); "
        f"catch e; disp(e.message); end; exit\""
    )

def main():
    # Parameter grids
    N_values        = [400]
    K_values        = [100]
    D_values        = [1]       # Average Diameter
    G_values        = [1.4]     # Ratio of large to small particles
    M_values        = [1]       # Mass of particles
    P_target_values = [.1, .01, .001]
    seed_values     = [1, 2, 3]
    x_mult_values   = [1]
    y_mult_values   = [1]
    z_mult_values   = [0]
    plotit          = False
    calc_eig        = False
    hertzian        = True
    save_path       = "./data/packings/hertz/"

    output_file = "./commandsPack.txt"

    with open(output_file, "w") as file:
        for N in N_values:
            for K in K_values:
                for D in D_values:
                    for G in G_values:
                        for M in M_values:
                            for P_target in P_target_values:
                                for seed in seed_values:
                                    for x_mult in x_mult_values:
                                        for y_mult in y_mult_values:
                                            for z_mult in z_mult_values:
                                                command = generate_matlab_command(
                                                        N, K, D, G, M,
                                                        P_target, seed,
                                                        plotit, x_mult, y_mult, z_mult,
                                                        calc_eig, save_path, hertzian
                                                        )
                                                file.write(command + "\n")

    print(f"Commands written to: {output_file}")

if __name__ == "__main__":
    main()

