# /Users/adams/anaconda3/envs/ann_solo/bin/python3

psms = pd.read_csv()

def get_mass_groups(psms, tol_mass, tol_mode, min_group_size=None):
    psms_remaining = psms.sort_values('search_engine_score[1]',
                                      ascending=False)
    psms_remaining['mass_diff'] = ((psms_remaining['exp_mass_to_charge'] -
                                    psms_remaining['calc_mass_to_charge']) *
                                   psms_remaining['charge'])

    # Start with the highest ranked SSM.
    mass_groups = []
    while psms_remaining.size > 0:
        # Find all remaining PSMs within the mass difference window.
        mass_diff = psms_remaining['mass_diff'].iat[0]
        if (tol_mass is None or tol_mode not in ('Da', 'ppm') or
                min_group_size is None):
            mask = np.full(len(psms_remaining), True, dtype=bool)
        elif tol_mode == 'Da':
            mask = (np.fabs(psms_remaining['mass_diff'] - mass_diff) <=
                    tol_mass)
        elif tol_mode == 'ppm':
            mask = (np.fabs(psms_remaining['mass_diff'] - mass_diffs) /
                    psms_remaining['exp_mass_to_charge'] * 10 ** 6
                    <= tol_mass)
        mass_groups.append(psms_remaining[mask])
        # Exclude the selected PSMs from further selections.
        psms_remaining = psms_remaining[~mask]

    mass_group_stats = []
    for mass_group in mass_groups:
        mass_group_stats.append((mass_group['mass_diff'].median(),
                                 mass_group['mass_diff'].mean(),
                                 len(mass_group)))
    mass_group_stats = pd.DataFrame.from_records(
        mass_group_stats, columns=['mass_diff_median', 'mass_diff_mean',
                                   'num_psms'])
    return mass_group_stats