import yaml
import argparse
import numpy as np
import multiprocessing as mp
from datetime import datetime

from NonlinearDissipativeSystems.src import Constants, RateCalc


if __name__ == '__main__':
    # Set up argument parser with all possible input parameters and defaults
    parser = argparse.ArgumentParser(description='Compute rates for a linear bath',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('-t', '--temperature', type=float, default=300.0, help='Temperature of the simulation in Kelvin')
    parser.add_argument('--n_bathmodes', type=int, default=9, help='Number of bathmodes')
    parser.add_argument('-g', '--gamma_factor', type=float, required=True, help="Friction value in units of omega_b")
    parser.add_argument('-n', '--n_beads', type=int, default=2, required=True, help="Number of beads for RPMD simulation")
    parser.add_argument('-l', '--linear', action='store_true', help="Flag to use linear bath coupling")
    parser.add_argument('-c', '--config', type=str, default='config.yaml', help="Path to config file")
    args = parser.parse_args()

    # Load config yaml file
    with open(args.config, 'r') as file:
        config = yaml.safe_load(file)

    n_samp_kappa = config['n_samp_kappa']
    n_samp_fe = config['n_samp_fe']
    n_samp_rd = config['n_samp_rd']
    n_equil = config['n_equil']
    n_evol_kappa = config['n_evol_kappa']
    n_evol_fe = config['n_evol_fe']
    n_points_fe = config['n_points_fe']
    n_repeats = config['n_repeats']
    save_path = config['save_path']

    ratecalc = RateCalc(
        T=args.temperature, 
        n_bathmodes=args.n_bathmodes, 
        n_samp_kappa=n_samp_kappa, 
        n_samp_fe=n_samp_fe,
        n_samp_rd=n_samp_rd,
        n_equil=n_equil, 
        n_evol_kappa=n_evol_kappa, 
        n_evol_fe=n_evol_fe,
        linear=args.linear
    )

    samples = np.linspace(1, n_repeats, n_repeats)
    linear_save_flag = 'l' if args.linear else 'nl'

    if args.n_beads == 1:
        # Checking if simulation is classical (i.e. 1 RP bead)
        #  If so, run the corresponding classical functions with Python multiprocessing

        pool = mp.Pool(mp.cpu_count())
        results_kappa_cl = pool.starmap(ratecalc.classical_transmission,
                                     [(args.gamma_factor, entry) for entry in samples])
        pool.close()

        # Label the output files with the current date/time
        now = datetime.now()
        dt_string = now.strftime("%d-%m-%Y_%H-%M-%S")
        np.savetxt(f'{save_path}/{linear_save_flag}_kappa_{dt_string}.csv', results_kappa_cl, delimiter=',')

    else:
        # Run RPMD rate calculation with Python multiprocessing

        pool = mp.Pool(mp.cpu_count())
        results_kappa = pool.starmap(ratecalc.rpmd_transmission,
                                     [(args.gamma_factor, args.n_beads, entry) for entry in samples])
        results_fe = pool.starmap(ratecalc.rpmd_mf_array,
                                  [(args.gamma_factor, args.n_beads, n_points_fe, entry) for entry in samples])
        results_rd = pool.starmap(ratecalc.reactant_distribution,
                                  [(args.gamma_factor, args.n_beads, entry) for entry in samples])

        pool.close()

        # Label the output files with the current date/time
        now = datetime.now()
        dt_string = now.strftime("%d-%m-%Y_%H-%M-%S")
        np.savetxt(f'{args.path_string}/{linear_save_flag}_kappa_gamma_{args.gamma_factor}.csv', results_kappa, delimiter=',')
        np.savetxt(f'{args.path_string}/{linear_save_flag}_mf_gamma_{args.gamma_factor}.csv', results_fe, delimiter=',')
        np.savetxt(f'{args.path_string}/{linear_save_flag}_rdist_gamma_{args.gamma_factor}.csv', results_rd, delimiter=',')