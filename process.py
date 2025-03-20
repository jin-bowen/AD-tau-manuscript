import pandas as pd
import re

def main():

	indel_mat = pd.read_csv('All_indel.tab', index_col=0, header=0)
	sig_mat = pd.read_csv('test.csv', header=0)
	indel_mat = indel_mat.T
	twobp_cols = [ bool(re.match('2bp_deletion_*', m)) for m in indel_mat.columns]
	twobp_cols = indel_mat.columns[twobp_cols]
	twobp_indel_mat = indel_mat[twobp_cols].sum(axis=1).to_frame('2bp_del')

	all_del_cols = [ bool(re.match(r'.*deletion*.', c)) for c in indel_mat.columns]	
	all_del_cols = indel_mat.columns[all_del_cols]
	all_del_indel_mat = indel_mat[all_del_cols].sum(axis=1).to_frame('all_del')
	meta_input = 'metafiles/all_meta.tab'
	meta = pd.read_csv(meta_input, header=0)

	tau_spl = meta.loc[meta['group']=='Tau','sample']
	notau_spl = meta.loc[meta['group']=='noTau', 'sample']
	sig_mat_tau = sig_mat.loc[sig_mat['sample'].isin(tau_spl)]	
	sig_mat_notau = sig_mat.loc[sig_mat['sample'].isin(notau_spl)]	
	
	all_tau = pd.merge(all_del_indel_mat.loc[tau_spl], twobp_indel_mat.loc[tau_spl], \
		left_index=True, right_index=True)
	all_tau = pd.merge(all_tau, sig_mat_tau, left_index=True, right_on='sample')
	all_notau = pd.merge(all_del_indel_mat.loc[notau_spl], twobp_indel_mat.loc[notau_spl], \
		left_index=True, right_index=True)
	all_notau = pd.merge(all_notau, sig_mat_notau, left_index=True, right_on='sample')
	
	all_tau[['sample','donor','burden','2bp_del','all_del','ID4']].to_csv('Tau_neuron.tab', index=False)
	all_notau[['sample','donor','burden','2bp_del','all_del','ID4']].to_csv('noTau_neuron.tab', index=False)


if __name__ == '__main__':
	main()

