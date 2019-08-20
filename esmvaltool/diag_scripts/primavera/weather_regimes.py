import os
import logging
import pickle

import matplotlib.pyplot as plt

import numpy as np
from scipy import stats

import iris
import iris.cube
import iris.analysis
import iris.util
import iris.coord_categorisation

import esmvaltool.diag_scripts.shared
import esmvaltool.diag_scripts.shared.names as n
from esmvaltool.diag_scripts.shared import group_metadata

from esmvaltool.diag_scripts.primavera.WRtool.CLUS_tool.WRtool.EOFtool import eof_computation
import esmvaltool.diag_scripts.primavera.WRtool.CLUS_tool.ctool as ctool

logger = logging.getLogger(os.path.basename(__file__))


class WeatherRegimes(object):
    def __init__(self, config):
        self.cfg = config
        self.filenames = esmvaltool.diag_scripts.shared.Datasets(self.cfg)
        self.variance_explained = self.cfg.get('variance_explained', None)
        self.pcs_selected = self.cfg.get('pcs', 4)

    def compute(self):
        data = group_metadata(self.cfg['input_data'].values(), 'alias')
        for alias in data:
            logger.info('Processing %s', alias)
            zg = iris.load_cube(data[alias][0]['filename'])
            logger.debug(zg)
            if not zg.coords('season'):
                iris.coord_categorisation.add_season(zg, 'time')
            for season in set(zg.coord('season').points):
                zg_season = zg.extract(iris.Constraint(season=season))
                logger.debug(zg_season)
                pc_unscaled = self._compute_eofs_and_pcs(alias, zg_season)
                self._clustering(alias, zg, pc_unscaled)

    def _compute_eofs_and_pcs(self, alias, zg):
        logger.debug('Compute and save EOFs and PCs')
        solver, pcs_scal1, eofs_scal2, pcs_unscal0, eofs_unscal0, varfrac = eof_computation(
            zg.data[np.newaxis, ...], zg.coord('latitude').points, zg.coord('longitude').points)

        # neof=0   # EOF to plot (neof starts from zero!)
        # tit='{0} {1} {2} {3} {4} {5}'.format(varname,model,kind,res,season,area)  # field decomposed with the EOF analysis
        #figPC_scal1, figEOF_scal2=eof_plots(neof,pcs_scal1, eofs_scal2,var,varunits,lat,lon,tit,numens)
        # plot the PCs
        # namef='{0}PCs{1}_{2}.eps'.format(OUTfig,neof+1,name_outputs)
        # figPC.savefig(namef)#bbox_inches='tight')
        # plot the EOFs
        # namef='{0}EOFs{1}_{2}.eps'.format(OUTfig,neof+1,name_outputs)
        # figEOF.savefig(namef)#bbox_inches='tight')
        #print('PCs and EOFs eps figure are saved in {0}'.format(OUTfig))
        # print('____________________________________________________________________________________________________________________')

        acc = np.cumsum(varfrac*100)
        if self.variance_explained:
            # Find how many PCs explain a certain percentage of variance
            # (find the mode relative to the percentage closest to perc, but bigger than perc)
            self.pcs_selected = min(
                enumerate(acc), key=lambda x: x[1] <= self.variance_explained)[0]+1
            logger.info(
                'The number of PCs that explain the percentage closest to'
                '{0}% of variance (but greater than {0}%) is {1}'.format(
                    self.variance_explained, self.pcs_selected
                ))
            exctperc = min(
                enumerate(acc), key=lambda x: x[1] <= self.variance_explained)[1]

        if self.pcs_selected:
            exctperc = acc[self.pcs_selected-1]
        logger.debug('(the first %s PCs explain exactly the %.2f of variance)',
                     self.pcs_selected, exctperc)

        # save python object solver
        ofile = os.path.join(self.cfg[n.WORK_DIR],
                             'solver_{0}.p'.format(alias))
        pickle.dump(solver, open(ofile, 'wb'), protocol=2)

        pc_coord = iris.coords.DimCoord(
            range(1, self.pcs_selected + 1),
            var_name='num',
            units='1.0',
        )

        # save EOF unscaled
        ofile = os.path.join(self.cfg[n.WORK_DIR],
                             'eof_{0}.nc'.format(alias))
        eof_unscaled = iris.cube.Cube(
            eofs_unscal0[:self.pcs_selected],
            var_name='eof',
            units=zg.units,
            dim_coords_and_dims=(
                (pc_coord, 0),
                (zg.coord('latitude'), zg.coord_dims('latitude')),
                (zg.coord('longitude'), zg.coord_dims('longitude')),
            )
        )
        logger.debug(eof_unscaled)
        iris.save(eof_unscaled, ofile)
        #tit='EOFunscal n{0}'.format(neof+1)
        # ____________Plot the lon-lat map (Orthographic Projection)
        #fig0 = plot_ortho(eofs_unscal0[0], lat_area, lon_area, clat=50, clon=0, tit=tit)
        # fig0.show()
        # plt.show(block=True)

        ofile = os.path.join(self.cfg[n.WORK_DIR],
                             'eof_scaled{0}.nc'.format(alias))
        cube = iris.cube.Cube(
            eofs_scal2[:self.pcs_selected],
            var_name='eof_scaled',
            units=zg.units,
            dim_coords_and_dims=(
                (pc_coord, 0),
                (zg.coord('latitude'), zg.coord_dims('latitude')),
                (zg.coord('longitude'), zg.coord_dims('longitude')),
            )
        )
        logger.debug(cube)
        iris.save(cube, ofile)
        #tit='EOFscal2 n{0}'.format(neof+1)
        # ____________Plot the lon-lat map (Orthographic Projection)
        #fig2 = plot_ortho(eofs_scal2[npc], lat_area, lon_area, clat=50, clon=0, tit=tit)
        # fig2.show()
        # plt.show(block=True)
        # save PC unscaled
        # [time x PCi]
        ofile = os.path.join(self.cfg[n.WORK_DIR],
                             'pc_{0}.nc'.format(alias))
        pcs_unscaled = iris.cube.Cube(
            pcs_scal1[:, :self.pcs_selected],
            var_name='pc',
            units=1.0,
            dim_coords_and_dims=(
                (pc_coord, 1),
                (zg.coord('time'), 0),
            )
        )
        iris.save(pcs_unscaled, ofile)

        # save PC scaled: PCs are scaled to unit variance (divided by the square-root of their eigenvalue)
        # [time x PCi]
        ofile = os.path.join(self.cfg[n.WORK_DIR],
                             'pc_scaled_{0}.nc'.format(alias))
        pcs_scaled = iris.cube.Cube(
            pcs_scal1[:, :self.pcs_selected],
            var_name='pc_scaled',
            units=1.0,
            dim_coords_and_dims=(
                (pc_coord, 1),
                (zg.coord('time'), 0),
            )
        )
        iris.save(pcs_scaled, ofile)
        # save varfrac: the fractional variance represented by each EOF mode
        ofile = os.path.join(self.cfg[n.WORK_DIR],
                             'varfrac_{0}.txt'.format(alias))
        np.savetxt(ofile, varfrac, fmt='%1.10f')

        return pcs_unscaled

    def _clustering(self, alias, zg, pc_unscaled):
        # load PCunscal
        # k-means analysis using the subset of PCs
        # ______________________________________
        print('number of clusters: {0}'.format(numclus))

        npart = 100
        varopt = []
        indcl0 = []
        centr = []
        print('npart ={0}'.format(npart))
        for clusters in range(2, 7):
            nfcl, indcl1, centr1, varopt1, iseed = ctool.cluster_toolkit.clus_opt(
                clusters, npart, pc_unscaled.data)
            # indcl1 starts from 1, indcl0 starts from 0
            indcl0.append(np.subtract(indcl1, 1))
            centr.append(centr1)
            varopt.append(varopt1)
        indcl = indcl0[numclus-2]
        centr = centr[numclus-2]

        # save cluster index
        namef = '{0}indcl_{1}clus_{2}.txt'.format(
            OUTtxt, numclus, name_outputs)
        np.savetxt(namef, indcl, fmt='%d')

        # save cluster centroids
        namef = '{0}centr_{1}clus_{2}.txt'.format(
            OUTtxt, numclus, name_outputs)
        np.savetxt(namef, centr)

        # save cluster optimal variance ratio (this is needed for significance computation: clusters_sig.py)
        namef = '{0}varopt_2to6clus_{2}.txt'.format(
            OUTtxt, numclus, name_outputs)
        np.savetxt(namef, varopt, fmt='%1.10f')

        # Cluster ordering in decreasing frequency
        # _______________________
        centrORD, indclORD = cluster_orderingFREQ(indcl, centr, numclus)

        # save cluster index
        namef = '{0}indclORD_{1}clus_{2}.txt'.format(
            OUTtxt, numclus, name_outputs)
        np.savetxt(namef, indclORD, fmt='%d')

        # save cluster centroids
        namef = '{0}centrORD_{1}clus_{2}.txt'.format(
            OUTtxt, numclus, name_outputs)
        np.savetxt(namef, centrORD)

        # COMPUTE AND SAVE CLUSTERS PATTERNS
        # _______________________
        # compute cluster patterns
        cluspattORD = compute_clusterpatterns(numclus, var_ensList, indclORD)
        # save cluster patterns
        varsave = 'cluspattern'
        ofile = '{0}cluspatternORD_{1}clus_{2}.nc'.format(
            OUTnc, numclus, name_outputs)
        save_N_2Dfields(lat_area, lon_area, np.array(
            cluspattORD), varsave, var_units, ofile)

        print(
            'Cluster pattern netCDF variable (ordered by decreasing frequency of occurrence) is saved as\n{0}'.format(ofile))
        print('____________________________________________________________________________________________________________________')

        print('\n******************************************************************************')
        print('END {0}'.format(sys.argv[0]))
        print('*********************************************************************************')

        # PLOT AND SAVE CLUSTERS
        # _______________________
        # plot the cluster patterns, all in one panel
        # namef='{0}indclORD_{1}clus_{2}.txt'.format(OUTtxt,numclus,name_outputs)
        # indcl=np.loadtxt(namef)
        # if enstoselect!='no':
        #    tit='{0} {1} {2}{3} {4} {5} {6} ({7}-{8})'.format(varname,model,kind,enstoselect,res,season,area,syr,eyr)
        # else:
        #    tit='{0} {1} {2} {3} {4} {5} ({6}-{7})'.format(varname,model,kind,res,season,area,syr,eyr)
        # ax,cluspatt=computeplot_clusters(area,lon_area,lat_area,numclus,numpcs,var_ensList,indcl,tit)
        #
        # namef='{0}clus_patterns_{1}clus_{2}.eps'.format(OUTfig,numclus,name_outputs)
        # ax.figure.savefig(namef)#bbox_inches='tight')
        #print('Clusters eps figure for {0} weather regimes is saved as\n{1}'.format(area,namef))
        # print('____________________________________________________________________________________________________________________')
        #
        # save cluster patterns
        # varsave='cluspattern'
        # ofile='{0}cluspattern_{1}clus_{2}.nc'.format(OUTnc,numclus,name_outputs)
        # save_N_2Dfields(lat_area,lon_area,np.array(cluspatt),varsave,'m',ofile)
        #print('Cluster pattern netCDF variable for {0} weather regimes is saved as\n{1}'.format(area,ofile))
        # print('____________________________________________________________________________________________________________________')
        # quick plot
        #tit='cluster pattern 1'
        # ____________Plot the lon-lat map (Orthographic Projection)
        #fig = plot_ortho(cluspatt[0], lat_area, lon_area, clat=50, clon=0, tit=tit)
        # fig.show()
        # plt.show(block=True)


def main():
    with esmvaltool.diag_scripts.shared.run_diagnostic() as config:
        WeatherRegimes(config).compute()


if __name__ == '__main__':
    main()
