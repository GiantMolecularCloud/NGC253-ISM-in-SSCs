#################################
# GAS IN SSCS: SPECTRAL FITTING #
#################################

# Get physical quantities from fits, such as integrated intensity, mass, linewidth, line shift, ...


###################################################################################################
# load data
###################################################################################################

execfile(os.path.join(basescriptdir, 'NGC253', 'paper_20a', '700.info.py'))
SSCs    = QTable.read(os.path.join(subprojectdir,'SSCs.fits'))
spectra = fnunpickle(os.path.join(specdir, 'spectra.pickle'))
bandfit = fnunpickle(os.path.join(specdir, 'band_fit_parameters.pickle'))
linedat = fnunpickle(os.path.join(specdir,'linedat.pickle'))

plotmols = [get_line(l) for l in ['CO 3-2 v=0', 'HCO+ 4-3', 'HCN 4-3', 'H13CN 4-3', 'HC15N 4-3', 'HC3N 39-38', 'CS 7-6', 'SO 8(8)-7(7) 3Sum_v=0', 'SO2 11(4,8)-11(3,9)']]


###################################################################################################
# generate postage stamps showing the distribution of gas around the SSCs
###################################################################################################

stamps = {str(SSC['no']): {} for SSC in SSCs}
for SSC in tqdm(SSCs):
    for plotmol in plotmols:
        try:

            # select correct frequency range
            mu     = linedat[str(SSC['no'])][plotmol['ID']]['mu']
            sigma  = linedat[str(SSC['no'])][plotmol['ID']]['sigma']

            # set region around SSC
            band   = 'LSB' if mu<350*u.GHz else 'USB'
            SSC_x  = SSC[band+' x']
            SSC_y  = SSC[band+' y']
            x0, x1 = [SSC_x-10,SSC_x+11]
            y0, y1 = [SSC_y-10,SSC_y+11]

            # collapse image
            bandim   = LSB if mu<350*u.GHz else USB
            specslab = bandim.spectral_slab(mu-2.5*sigma, mu+2.5*sigma)
            subcube  = specslab[:, x0:x1, y0:y1]
            mom0     = subcube.with_spectral_unit(u.km/u.s, velocity_convention='optical', rest_value=plotmol['restfreq']).moment0()

            stamps[str(SSC['no'])][plotmol['ID']] = {'moment 0': mom0}

        except:
            print("SSC "+str(SSC['no'])+": "+plotmol['ID']+" not fitted.")


###################################################################################################
# plot postage stamps showing the distribution of gas around the SSCs
###################################################################################################

# When plotting interactively, the image size are all wrong!
plt.ioff()

# set up image grid
nrows = len(SSCs)
ncols = len(plotmols)
top, bottom, left, right = [0.80,0.05,0.08,0.95]

fig,axes = plt.subplots(nrows=nrows, ncols=ncols, squeeze=True, sharex='col', sharey='row', figsize=(ncols-1.2,nrows))
fig.subplots_adjust(hspace=0, wspace=0, top=top, bottom=bottom, left=left, right=right)
# caxes = [fig.add_axes([left+i_mol*(right-left)/ncols, top, (right-left)/ncols, 0.02]) for i_mol,_ in enumerate(plotmols)]

# run loop to plot stamps
for i_SSC,SSC in tqdm(enumerate(SSCs)):
    for i_mol,(plotmol,scale) in enumerate(zip(plotmols,[7000,1500,1500,250,250,500,1500,500,500])):
        i_im = i_SSC*len(SSCs)+i_mol+1
        ax = axes[i_SSC][i_mol]

        try:
            mom0 = stamps[str(SSC['no'])][plotmol['ID']]['moment 0']
            # show map
            im = ax.imshow(mom0.value, origin='lower', cmap='binary', vmin=-0.1*scale, vmax=scale)

            # show SSC contour
            cdx = u.Quantity(str(mom0.header['cdelt1'])+mom0.header['cunit1'])
            cdy = u.Quantity(str(mom0.header['cdelt2'])+mom0.header['cunit2'])
            x = int(np.floor(mom0.shape[0]/2))
            y = int(np.floor(mom0.shape[1]/2))
            FWHM_deg = (np.sin( (SSC['FWHM']/distance).to(u.dimensionless_unscaled).value )*u.radian).to(u.degree)
            w = (FWHM_deg/cdx).value
            h = (FWHM_deg/cdy).value
            a = 0.
            ax.add_patch(mpl.patches.Ellipse(xy=[x,y], width=w, height=h, angle=a, facecolor='none', edgecolor='r',  linewidth=1.))
        except:
            continue
        finally:
            # panel formatting
            ax.set_aspect('equal')
            major_ticks = [2,6,10,14,18]
            ax.set_xticks(major_ticks)
            ax.set_yticks(major_ticks)
            ax.tick_params(labelleft=False, labelright=False, labeltop=False, labelbottom=False, left=True, right=True, top=True, bottom=True)
            # if i_SSC==0:
            #     fig.colorbar(im, cax=caxes[i_mol], orientation='horizontal')
            #     caxes[i_mol].xaxis.set_ticks_position('top')

# show labels in bottom left panel
# axes[-1][0].xaxis.set_visible(True)
axes[-1][-1].tick_params(labelbottom=True, labelright=True, bottom=True, right=True)
axes[-1][-1].yaxis.set_label_position('right')
axes[-1][-1].set_xticklabels(['-0.4','-0.2','0.0','0.2','0.4'], fontsize=8, rotation=90)
axes[-1][-1].set_yticklabels(['-0.4','-0.2','0.0','0.2','0.4'], fontsize=8)
axes[-1][-1].set_xlabel(r'offset [$^{\prime\prime}$]', fontsize=8)
axes[-1][-1].set_ylabel(r'offset [$^{\prime\prime}$]', fontsize=8)

# label SSCs
for i_SSC,SSC in tqdm(enumerate(SSCs)):
    ax = axes[i_SSC][0]
    ax.text(-0.1, 0.5, 'SSC '+str(SSC['no']), color='k', fontsize=12, weight='bold', transform=ax.transAxes, ha='right', va='center')

# label lines
for i_mol,plotmol in tqdm(enumerate(plotmols)):
    ax = axes[0][i_mol]
    ax.text(0.5, 1.1, line_tex(plotmol), color='k', fontsize=12, weight='bold', transform=ax.transAxes, ha='center', va='bottom', rotation=90)

# fig.tight_layout()
fig.savefig(os.path.join(plotdir, 'stamps.pdf'), dpi=300, bbox_inches='tight')


###################################################################################################
#
###################################################################################################
