import numpy as np

def generate_matlab_command(N, K, D, G, M, P_target, seed, plotit, x_mult, y_mult, calc_eig, save_path):
    plotit_str = "true" if plotit else "false"
    calc_eig_str = "true" if calc_eig else "false"
    return (
        f"matlab -nodisplay -nosplash -r \"addpath('./src/matlab_functions/'); try; "
        f"pack({N}, {K}, {D}, {G}, {M}, {P_target}, {seed}, {plotit_str}, {x_mult}, {y_mult}, 0, {calc_eig_str}, '{save_path}'); catch e; disp(e.message); end; exit\""
    )

def main():
    # Define different values for each variable
    N_values = [400]
    K_values = [100]
    D_values = [1]           # Average Diameter
    G_values = [1.4]         # Ratio of large to small particles
    M_values = [1]           # Mass of particles
    P_target_values = [.0003162, .003162, .03162]
    #P_target_values = np.logspace(-4, -1,10)
    seed_values = [1, 2, 3, 4, 5]
    plotit = False
    x_mult_values = [400]
    y_mult_values = [1]
    calc_eig = False
    save_path = './in/2d/10by10/'
    output_file = "./commandsPack2d.txt"

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
                                            command = generate_matlab_command(
                                                N, K, D, G, M, P_target, seed,
                                                plotit, x_mult, y_mult, calc_eig,
                                                save_path
                                            )
                                            file.write(command + "\n")
    print(f"Commands written to: {output_file}")

if __name__ == "__main__":
    main()


