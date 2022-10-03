import pandas as pd
import numpy as np


class CAMML_class:
	def __init__(self, filepath, mdlpath):
		sheet_names = ['Baseline', 'Drilling','Hydraulic Fracturing', 'Mill Out', 'FlowbackProduction']
		print('CAMML upload started')

		# import data
		# VOC data taken every hour, contains flags for instrument mishaps -- see flag file in data folder
		self.imp_voc_data = pd.concat([pd.read_excel(filepath + '\Extraction_Livingston_VOC_Data.xlsx', sheet_name=s,
													 engine='openpyxl').assign(sheet_name=s) for s in sheet_names])
		self.imp_voc_data.reset_index(inplace=True, drop=True)

		# Other data includes O3, PM, NOx, etc taken every minute
		self.imp_other_data = pd.concat([pd.read_excel(filepath + '\Extraction_Livingston_Other_Data.xlsx',
													   sheet_name=s, engine='openpyxl').assign(sheet_name=s) for s in sheet_names])
		self.imp_other_data.reset_index(inplace=True, drop=True)

		# takes cos of wind direction to deal with circular data
		self.imp_other_data['cosine'] = np.cos(np.radians(self.imp_other_data['WindDirection_[degrees]']))
		self.imp_other_data['sin'] = np.sin(np.radians(self.imp_other_data['WindDirection_[degrees]']))

		# CAMML method detection limit
		self.mdl = pd.read_excel(mdlpath + '\CAMML_Teir1_202008_MDLs.xlsx', index_col='Compound', engine='openpyxl')

		# create new df for data vs flags
		self.voc_data = pd.DataFrame()
		self.voc_data_flags = pd.DataFrame()
		self.other_data = pd.DataFrame()
		self.other_data_flags = pd.DataFrame()
		self.unc = pd.DataFrame(index=range(3780))

		self.merged_data = pd.DataFrame()

		# good to knows
		self.voccols = []
		self.othercols = []
		self.allcols = []



	def camml_clean(self):
		#  splitting voc into data and flags
		self.imp_voc_data['midtime'] = [
			self.imp_voc_data.loc[time, 'StartTime_[LST]'] + ((self.imp_voc_data.loc[time, 'EndTime_[LST]'] - self.imp_voc_data.loc[time, 'StartTime_[LST]']) / 2)
			for time in self.imp_voc_data.index]

		self.voc_data_flags = pd.DataFrame(self.imp_voc_data.loc[:, self.imp_voc_data.columns.str.contains('Flags')])  # creating flag df
		self.voc_data = pd.DataFrame(self.imp_voc_data.loc[:, ~self.imp_voc_data.columns.str.contains('Flags')])  # creating concentration df
		self.voc_data.set_index('midtime', inplace=True)

		self.voc_data_flags.fillna(9, inplace=True)  # Data and flags filled
		self.voc_data_flags['midtime'] = self.voc_data.index
		self.voc_data_flags.set_index('midtime', inplace=True)

		#  splitting other into data and flags
		self.imp_other_data['midtime'] = [self.imp_other_data.loc[time, 'StartTime_[LST]'] +
										  ((self.imp_other_data.loc[time, 'EndTime_[LST]'] -
											self.imp_other_data.loc[time, 'StartTime_[LST]']) / 2) for time in
										  self.imp_other_data.index]
		self.other_data_flags = pd.DataFrame(
			self.imp_other_data.loc[:, self.imp_other_data.columns.str.contains('Flags')])  # creating flag df
		self.other_data = pd.DataFrame(
			self.imp_other_data.loc[:, ~self.imp_other_data.columns.str.contains('Flags')])  # creating concentration df
		self.other_data.set_index('midtime', inplace=True)

		self.other_data_flags.fillna(9, inplace=True)  # Data and flags filled
		self.other_data_flags['midtime'] = self.other_data.index
		self.other_data_flags.set_index('midtime', inplace=True)
		print('finished camml_clean')

	def voc_flags(self):
		self.voc_data.drop(['2-methyl-pentane_[ppbV]', '3-methyl-pentane_[ppbV]',
							'a-pinene_[ppbV]'], axis=1,inplace=True)  # These are columns that contain too many 2,3 flags
		self.voc_data_flags.drop(['2-methyl-pentane_Flags', '3-methyl-pentane_Flags',
							'a-pinene_Flags'], axis=1,inplace=True)
		self.voc_data.drop(['p-ethyl-toluene_[ppbV]', '1,2,4-trimethyl-benzene_[ppbV]',
							'p-diethyl-benzene_[ppbV]'], axis=1, inplace=True)  # these are columns that have greater than 85% of data below LOD (raised 1 or 4 flag)
		self.voc_data_flags.drop(['p-ethyl-toluene_Flags', '1,2,4-trimethyl-benzene_Flags',
							'p-diethyl-benzene_Flags'], axis=1, inplace=True)

		# now I need to replace values in remaining dataset
		for mdl, col in zip(self.mdl, self.mdl.index):
			try:
				self.voc_data[col].where(self.voc_data[col] > mdl, mdl/2, inplace=True)  # this takes care of 1, 4 flags (not actually b/c 0.3 raised flags but for my purposes)
			except KeyError:
				continue

		# Now I need to take care of 2,4 flags
			flag23 = ((self.voc_data_flags == 2)|(self.voc_data_flags == 3))[col+'_Flags'].values  # locations of 2,3 flags in flag dataset
			dropme = ((self.voc_data_flags == 2)|(self.voc_data_flags == 3))[col+'_Flags'].index
			# replacing flag23 values with median of other values to take care of 2,3 flags
			self.voc_data.loc[flag23, col+['_[ppbV]']] = self.voc_data.loc[[i for i in self.voc_data.index if i not in dropme],col+['_[ppbV]']].median()



	def avg_other(self):
	#  Averaging other data to fit into voc data time ranges
		for idx in self.voc_data.index:  # Taking other data, finding times between and adding idx column
			self.other_data.loc[(self.voc_data.loc[idx, 'StartTime_[LST]'] < self.other_data.index) & (
						self.voc_data.loc[idx, 'EndTime_[LST]'] > self.other_data.index), 'between_time'] = idx


		self.other_data.columns = [i[0] for i in self.other_data.columns.str.split('_')]
		self.other_data.rename(columns={'PM2.5(LTP)':'PM2.5', 'PM10(LTP)':'PM10', 'between':'between_time'}, inplace=True)

		self.othercols = [i for i in self.other_data.columns if any(k in i for k in ['PM', 'CH', 'NO', 'O3','sin','cosine'])]


		for num, col in enumerate(self.othercols):
			k = self.other_data[['between_time', col]].groupby('between_time').agg('mean')
			self.merged_data[col] = k[col]


		# merging data
		self.merged_data = pd.merge_asof(self.voc_data, self.merged_data, left_index=True, right_index=True, direction='nearest')  # changed right on to right index
		self.merged_data['PM2.510'] = self.merged_data['PM10'] - self.merged_data['PM2.5']
		self.merged_data['WindDirection'] = np.arctan2(self.merged_data['sin'],self.merged_data['cosine']) * (180/np.pi)


	def unc_func(self):
		print('starting unc_func')
		# rename merged data columns
		self.merged_data.columns = [i[0] for i in
									self.merged_data.columns.str.split('_')]  # getting rid of everything after _


		# column masks
		self.voccols = [i for i in self.voc_data.select_dtypes(include=float).columns if
						i != 'a-pinene_[ppbV]']  # all the vocs with data
		self.voccols = [i.replace('_[ppbV]', '') for i in self.voccols]  # is this actually gonna work...?

		self.othercols = [i for i in self.merged_data.columns if any(k in i for k in ['PM', 'CH', 'NO', 'O3'])]

		self.allcols = self.voccols + self.othercols
		# doubling unc for OTHER cols
		othermdl = 0.03  # self.merged_data[self.othercols].quantile(0.25) * 0.05  # mdl for other compounds
		self.unc.index = self.merged_data.index.copy(deep=True)
		# mean values for unc, base unc based on 5% of mean concentration or 0.05 whichever is larger
		# Polissar analytical unc + lod/3 as determined data unc
		for col in self.allcols:  # create df where unc == [] * 5% + lod/3
			self.merged_data.loc[:,col] = self.merged_data.loc[:,col].fillna(-999)
			self.unc.loc[:,col] = self.merged_data.loc[:,col].copy(deep=True) * 0.05 # unc is 5%

			if col in self.voccols:
				k = lambda x: x + (self.mdl.loc[col].item() / 3)  # LOD is given by self.mdl
				self.unc[col] = self.unc[col].apply(k).copy(deep=True)  # unc is now unc * conc + lod/3

			# if data below lod,data == lod/2 and unc == 5/6 lod.
			# determined values have unc == unc * [] + lod/3
			# where lod == self.mdl and unc == 5%
			# Now within the VOC cols I need to find [] that are below the lod

				below_mdl_mask = self.merged_data[col] < self.mdl.loc[col].item()
				idx = self.merged_data[below_mdl_mask].index.copy(deep=True)
				self.unc.loc[idx, col] = (5/6) * np.copy(self.mdl.loc[col].item())


			elif col in self.othercols:
				q = lambda u: u + othermdl / 3  # LOD is 0.03
				self.unc[col] = self.unc[col].apply(q).copy(deep=True)  # unc is now unc * conc + lod/3

				below_mdl_mask = self.merged_data[col] < othermdl
				idx = self.merged_data[below_mdl_mask].index.copy(deep=True)
				self.unc.loc[idx, col] = (5 / 6) * othermdl

		# if data below lod,data == lod/2 and unc == 5/6 lod.
		# determined values have unc == unc * [] + lod/3
		# where lod == self.mdl and unc == 5%
		# Now I need to replace the data within the data to be 1/2 mdl
		[self.merged_data[col].where(self.merged_data[col] > self.mdl.loc[col].item(), self.mdl.loc[col].item() / 2,
									 inplace=True) for col in self.voccols]  # only replaced VOC data
		[self.merged_data[col].where(self.merged_data[col] > othermdl, othermdl / 2,
									 inplace=True) for col in self.othercols]



	def flow(self):
		print('running flow')
		self.camml_clean()
		self.voc_flags()
		print('finished VOC flags')
		self.avg_other()
		print('avgeraging other finished')
		self.unc_func()
		print('done!')















