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


def simulate_effect_gamma(
    selected_variants, n_var, effect_shape=0.4, effect_scale=1.66
):
    # pull the effect size from a gamma distribution
    effects = np.random.default_rng().gamma(effect_shape, effect_scale, n_var)

    directions = np.random.choice([-1, 1], n_var)
    effects = effects * directions

    # Create a DataFrame from the selected variants array and effects
    result_df = pd.DataFrame({"Marker": selected_variants, "EFFECT": effects})

    return result_df


def simulate_effect_uniform(selected_variants, n_var, low_end=0.1, high_end=0.5):
    # pull the effect size from a uniform distribution
    effects = np.random.uniform(low_end, high_end, n_var)

    directions = np.random.choice([-1, 1], n_var)
    effects = effects * directions

    # Create a DataFrame from the selected variants array and effects
    result_df = pd.DataFrame({"Marker": selected_variants, "EFFECT": effects})

    return result_df


if __name__ == "__main__":
    # Get command line arguments
    strain_set_variant_file = sys.argv[1]
    n_var = int(sys.argv[2])

    if sys.argv[3] == "gamma":
        effect_type = "gamma"
    elif "-" in sys.argv[3]:
        effect_type = "uniform"
        effect_range = sys.argv[3].split("-")
        low_end = float(effect_range[0])
        high_end = float(effect_range[1])
    else:
        print(
            f"Error: Invalid effect specification '{sys.argv[3]}'. Use 'gamma' or 'low_end-high_end'"
        )
        sys.exit(1)
    # Read in annotated strain_set variants
    strain_var = load_strain_set_variants(strain_set_variant_file)

    # Select Causal variants for orthogroups
    causal_vars = select_variants(strain_var, n_var)
    print("Selected Causal Variants")

    if effect_type == "gamma":
        print("Simulating effects using gamma distribution")
        causal_vars_effects = simulate_effect_gamma(causal_vars, n_var)
    else:
        print(f"Simulating effects using uniform distribution: {low_end}-{high_end}")
        causal_vars_effects = simulate_effect_uniform(
            causal_vars, n_var, low_end=low_end, high_end=high_end
        )
    print("Simulated effects")

    # Write output for trait simulations - just id and effect
    causal_vars_effects.to_csv(
        "causal_vars.txt",
        sep=" ",
        index=False,
        header=False,
    )
