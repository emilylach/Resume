import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt


class PMFAnalysis:

    def __init__(self, path, fileprefix, factor_num=6, sam_len=1131, compound_len=34):
        self.sam_len = sam_len
        self.factor_num = factor_num
        self.compound_len = compound_len
        self.contributions = pd.read_csv(path.strip() + fileprefix + '_contributions.csv', header=1135)
        self.profiles = pd.read_csv(path.strip() + fileprefix + '_profiles.csv')
        self.conc_prof = pd.DataFrame()
        self.perc_species_prof = pd.DataFrame()
        self.perc_total_prof = pd.DataFrame()


        self.withsite = pd.DataFrame()
        self.construction_data = pd.read_excel(r'C:\Users\emlach\Documents\Masters Thesis\Broomfield\Excel Data\CCOB_construction.xlsx', engine='openpyxl')
        self.construction_group = self.construction_data.groupby('Pad')
        self.concentrations_by_compound = np.zeros((self.factor_num, self.sam_len, self.compound_len))

    def contributions_func(self):
        """ Contributions is contributions to Factors. Makes the timeline series from PMF """
        self.contributions.reset_index(inplace=True)
        self.contributions.drop('Unnamed: 0', axis=1, inplace=True)
        self.contributions.drop(0, inplace=True)
        self.contributions['middate'] = [dt.datetime.strptime(self.contributions.loc[k, 'Unnamed: 1'], '%m/%d/%y %H:%M')
                                         for k in self.contributions.index]
        self.contributions.drop('Unnamed: 1', axis=1, inplace=True)
        self.contributions.set_index('middate', inplace=True)
        # self.contributions.columns = ['Factor %s' % i for i in range(1, self.factor_num + 1)]
        self.contributions = self.contributions.astype(float)

    def profiles_func(self):
        self.profiles.reset_index(inplace=True)
        self.profiles.drop(self.profiles.iloc[:,[0,1]].columns, axis=1, inplace=True)
        self.profiles.columns = ['compounds']+ [f'Factor {i}' for i in range(1,7)]
        self.profiles.drop([0,1], inplace=True)
        self.profiles.reset_index(inplace=True, drop=True)
        self.profiles.drop(self.profiles.loc[34:, :].index)

        self.conc_prof = self.profiles.drop(self.profiles.loc[self.compound_len:, :].index)
        self.conc_prof.set_index('compounds', inplace=True)
        self.conc_prof = self.conc_prof.astype(float)

        self.perc_species_prof = self.profiles.loc[34 + 4: 2 * 34 + 3, :].copy(deep=True)
        self.perc_species_prof.set_index('compounds', inplace=True)
        self.perc_species_prof = self.perc_species_prof.astype(float)
        self.perc_species_prof = self.perc_species_prof/100


        self.perc_total_prof = self.profiles.loc[2*self.compound_len + 8:, :].copy(deep=True)
        self.perc_total_prof.set_index('compounds', inplace=True)
        self.perc_total_prof = self.perc_total_prof.astype(float)
        self.perc_total_prof = self.perc_total_prof/100

    def inratio(self):
        return dict(zip(self.conc_prof.columns, [self.conc_prof.loc['i_pentane', fac]/self.conc_prof.loc['n_pentane', fac] for fac in self.conc_prof.columns]))

    def tbratio(self):
        return dict(zip(self.conc_prof.columns, [self.conc_prof.loc['toluene', fac]/self.conc_prof.loc['benzene', fac] for fac in self.conc_prof.columns]))

    def mergesite(self, df_wsite):
        """ Will merge on midtime. Will call srt_site from Clean Class """
        self.withsite = pd.merge_asof(self.contributions, df_wsite[['srt_site', 'midtime']].sort_values('midtime'), left_index=True, right_on='midtime')

    def grabsite(self, site_str):
        return self.withsite.loc[self.withsite['srt_site'].isin(site_str), :]

    def fac_plot(self, df, site, com_or_fac):
        fig, ax = plt.subplots()

        for i in self.construction_group.get_group('Livingston').index:
            ax.fill_betweenx(np.arange(0, np.amax(df[com_or_fac]) + 0.2, 0.01), self.construction_group
                             .get_group('Livingston')['Start Date'][i], self.construction_group.get_group('Livingston')
                            ['Stop Date'][i], alpha=0.5, label=self.construction_data['Activity'][i])

    def factor_by_compound(self):
        cols = []
        for fac in range(1, self.factor_num + 1, 1):
            sam_by_fac = self.contributions.loc[:, f'Factor {fac}'].values.reshape(self.sam_len, 1)
            fac_by_com = self.perc_total_prof[f'Factor {fac}'].T.values.reshape(1, self.compound_len)
            self.concentrations_by_compound[fac-1, :, :] = (sam_by_fac.dot(fac_by_com))  # sam_by_fac total concentration per compound over sample one
            for com in self.perc_species_prof.index:
                cols.append(f'{com} fac {fac}')





	# PMF 6 Factor analysis
	def imp_ccob6(self):
		self.contributions_func()
		self.profiles_func()
		self.mergesite()  # added in data before pmf

		self.fac_names = dict(zip([6, 3, 5, 1, 2, 4],
								  ['Light Alkanes', 'Ethyne', 'Background', 'Drilling', 'Complex Alkanes',
								   'Combustion']))



		def drop(df):
			dropidx = df[df['Species'].isnull()].index
			df.drop(dropidx, inplace=True)
			mask = df['Species'].str.contains('Factor')
			splitidx = list(df[mask].index)
			splitidx.append(df.iloc[-1].name + 1)
			return [df.loc[splitidx[i]:splitidx[i + 1] - 1].reset_index(drop=True).drop([0, 1]).set_index(
				'Species') for
					i in range(len(splitidx) - 1)]


		return self.ccob6