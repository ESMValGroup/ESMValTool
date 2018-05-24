"""Module for plotting hydrocycle metrics"""

import csv

GLOBAL_METRICS = ['Precipitation Global {}',
                  'Precipitation Land {}',
                  'Precipitation Ocean {}',
                  'Evaporation Global {}',
                  'Evaporation Land {}',
                  'Evaporation Ocean {}',
                  'Surface Runoff Land {}',
                  'Sub-surface Runoff Land {}',
                  'Total Runoff Land {}',
                  'Sensible heat flux Global {}',
                  'Sensible heat flux Land {}',
                  'Sensible heat flux Ocean {}',
                  'Latent heat flux Global {}',
                  'Latent heat flux Land {}',
                  'Latent heat flux Ocean {}',
                  'Total water vapour Global {}',
                  'Total water vapour Land {}',
                  'Total water vapour Ocean {}',
                  'Total cloud liquid water Global {}',
                  'Total cloud liquid water Land {}',
                  'Total cloud liquid water Ocean {}',
                  'Total cloud ice water Global {}',
                  'Total cloud ice water Land {}',
                  'Total cloud ice water Ocean {}']

REGIONS_METRICS = {'tropics': [
                      'Precipitation Global {}',
                      'Mean Precip: Tropics {}',
                      'Mean Evap: Tropics {}',
                      'Mean Water Vapour: Tropics {}',
                      'Mean Cloud liquid water: Tropics {}',
                      'Mean Cloud ice water: Tropics {}',
                      'Mean Sfc Latent heat: Tropical Ocean {}',
                      'Mean Sfc Sensible heat: Tropical Land {}',
                      'Total runoff: Tropical Land {}'],
                  'NHExtratropics': [
                      'Precipitation Global {}',
                      'Mean Precip: NH Extra-tropics {}',
                      'Mean Evap: NH Extra-tropics {}',
                      'Mean Water Vapour: NH Extra-tropics {}',
                      'Mean Cloud liquid water: NH Extra-tropics {}',
                      'Mean Cloud ice water: NH Extra-tropics {}',
                      'Mean Sfc Latent heat: NH Extra-tropical Ocean {}',
                      'Mean Sfc Sensible heat: NH Extra-tropical Land {}',
                      'Total runoff: NH Extra-tropical Land {}'],
                  'SHExtratropics': [
                      'Precipitation Global {}',
                      'Mean Precip: SH Extra-tropics {}',
                      'Mean Evap: SH Extra-tropics {}',
                      'Mean Water Vapour: SH Extra-tropics {}',
                      'Mean Cloud liquid water: SH Extra-tropics {}',
                      'Mean Cloud ice water: SH Extra-tropics {}',
                      'Mean Sfc Latent heat: SH Extra-tropical Ocean {}',
                      'Mean Sfc Sensible heat: SH Extra-tropical Land {}',
                      'Total runoff: SH Extra-tropical Land {}']}



def create_plot_selections():
    '''
    Routine to plot hydrocycle metrics in small groups
    '''
    seasons = ['ann', 'djf', 'mam', 'jja', 'son']

    # Metrics to match the old-style plots
    for season in seasons:
        metrics = [metric.format(season) for metric in GLOBAL_METRICS]
        write_csv('hydrocycle_global_' + season + '.csv', metrics)

    # Regions
    regions = ['tropics', 'NHExtratropics', 'SHExtratropics']
    for season in seasons:
        for region in regions:
            metrics = [metric.format(season) for metric in REGIONS_METRICS[region]]
            write_csv('hydrocycle_' + region + '_' + season + '.csv', metrics)

    # P minus E
    regions = ['Global', 'Land', 'Ocean']
    metrics = ['P minus E {} {}'.format(region, season) for season in seasons for region in regions]
    write_csv('hydrocycle_' + 'P_minus_E_' + 'all_seasons' + '.csv', metrics)


    # Interannual variations
    regions = ['Global', 'Tropics', 'NH Extra-tropics', 'SH Extra-tropics']
    metrics = ['St dev Precip: {} {}'.format(region, season) for region in regions for season in seasons]
    write_csv('hydrocycle_' + 'interannual_STDEV_' + 'all_seasons' + '.csv', metrics)


def write_csv(name, rows):
    with open(name, 'w') as fh:
        writer = csv.writer(fh)
        for row in rows:
            writer.writerow([row])


if __name__ == '__main__':
    create_plot_selections()
