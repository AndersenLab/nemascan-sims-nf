#!/usr/bin/env python3
import pandas as pd
import numpy as np
import sys


def load_strain_set_variants(strain_set_variant_file):
    """
    Load the strain set variants from a plink .bim file
    """

    # just read in the marker column from the .bim file
    strain_set_variants = pd.read_csv(
        strain_set_variant_file,
        sep="\t",
        header=None,
        usecols=[3],
        names=["Marker"],
    )

    return strain_set_variants


def select_variants(markers, n_var):
    """
    Select n_var variants from the list of markers
    """
    # Select n_var variants from the list of markers
    # Ensure markers is a 1D array by extracting the "Marker" column if it's a DataFrame
    if isinstance(markers, pd.DataFrame):
        markers = markers["Marker"].values

    selected_variants = np.random.choice(markers, n_var, replace=False)
    return selected_variants


# - [ ] allow dynamic effect shape & scale input
def simulate_og_effect_gamma(
    selected_variants, n_var, effect_shape=0.4, effect_scale=1.66
):
    # pull the effect size from a gamma distribution
    effects = np.random.default_rng().gamma(effect_shape, effect_scale, n_var)

    directions = np.random.choice([-1, 1], n_var)
    effects = effects * directions

    # Create a DataFrame from the selected variants array and effects
    result_df = pd.DataFrame({"Marker": selected_variants, "EFFECT": effects})

    return result_df


if __name__ == "__main__":
    # Define orthogroups from command line argument where the orthogroups are separated by commas

    # Get the list of variants in the strain sets - from the .bim file
    # strain_set_variant_file = sys.argv[1]
    strain_set_variant_file = (
        "/Users/ryanmckeown/Desktop/nemascan-sims-nf/data/plink/TO_SIMS.bim"
    )

    # get the number of varints to select from the command line argument
    # n_var = int(sys.argv[2])
    n_var = 100

    # Read in annotated strain_set variants
    strain_var = load_strain_set_variants(strain_set_variant_file)

    # Select Causal variants for orthogroups
    causal_vars = select_variants(strain_var, n_var)
    print("Selected Causal Variants")

    # Simulate effects for causal variants
    causal_vars_effects = simulate_og_effect_gamma(causal_vars, n_var)
    print("Simulated Effects for Causal Variants")

    # Write output for trait simulations - just id and effect
    causal_vars_effects.to_csv(
        "/Users/ryanmckeown/Desktop/nemascan-sims-nf/data/causal_og_vars.txt",
        sep=" ",
        index=False,
        header=False,
    )
