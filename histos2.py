from ROOT import TH1F


def histos():
	"""
	TH1F objetcs creation.
	- Convention of naming is: plot_kinematicOrTopologicalVariable_particle.
	- b1 an b2 are the leading b-jets and so on.
	Return: List of TH1F objetcs. 
	The last histo is plot_MET and the previous one is plot_PT_jets.
	"""
	#b-quark plots
	plot_PT_b1 = TH1F("PT_b1", "PT_b1", 100, 0.0, 1000.0)
	plot_ETA_b1 = TH1F("Eta_b1", "Eta_b1", 100, -5, 5) 
	plot_PHI_b1 = TH1F("Phi_b1", "Phi_b1", 100, -4, 4)  

	#tau plots
	plot_PT_tau1 = TH1F("PT_tau1", "PT_tau1", 100, 0.0, 1000.0) 
	plot_ETA_tau1 = TH1F("Eta_tau1", "Eta_tau1", 100, -5, 5)
	plot_PHI_tau1 = TH1F("Phi_tau1", "Phi_tau1", 100, -4, 4)
	plot_PT_tau2 = TH1F("PT_tau2", "PT_tau2", 100, 0.0, 1000.0) 
	plot_ETA_tau2 = TH1F("Eta_tau2", "Eta_tau2", 100, -5, 5) 
	plot_PHI_tau2 = TH1F("Phi_tau2", "Phi_tau2", 100, -4, 4)

	#jet plots
	plot_PT_j1 = TH1F("PT_jet1", "PT_jet1", 100, 0.0, 1000.0) 
	plot_ETA_j1 = TH1F("Eta_jet1", "Eta_jet1", 100, -5, 5)
	plot_PHI_j1 = TH1F("Phi_jet1", "Phi_jet1", 100,-4, 4)
	plot_PT_j2 = TH1F("PT_jet2", "PT_jet2", 100, 0.0, 1000.0) 
	plot_ETA_j2 = TH1F("Eta_jet2", "Eta_jet2", 100, -5, 5)
	plot_PHI_j2 = TH1F("Phi_jet2", "Phi_jet2", 100,-4, 4)

	#Delta R
	plot_DeltaR_tau1tau2 = TH1F("DeltaR_tau1tau2", "DeltaR_tau1tau2", 100, 0, 10)
	plot_DeltaR_j1j2 = TH1F("DeltaR_j1j2", "DeltaR_j1j2", 100, 0,10)

	#Delta Phi
	plot_DeltaPhi_tau1tau2 = TH1F("DeltaPhi_tau1tau2", "DeltaPhi_tau1tau2", 100, -4, 4)
	plot_DeltaPhi_j1j2 = TH1F("DeltaPhi_j1j2", "DeltaPhi_j1j2", 100, -4, 4)
	
	#Delta Pt
	plot_DeltaPt_tau1tau2 = TH1F("DeltaPt_tau1tau2", "DeltaPt_tau1tau2", 100, 0, 10)
	plot_DeltaPt_j1j2 = TH1F("DeltaPt_j1j2", "DeltaPt_j1j2", 100, 0,10)
	
	#Reconstructed masses (j,j)
	plot_M_Zp = TH1F("M_Zp", "M_Zp", 100, 0.0, 1500.0)
	plot_M_j1j2 = TH1F("M_j1j2", "M_j1j2", 100, 0.0, 1000.0)
	
	#Reconstructed masses (j,j,b1)
	plot_M_j1j2b1 = TH1F("M_j1j2b1", "M_j1j2b1", 100, 0.0, 1000.0)
	
	#StMet
	plot_StMet = TH1F("StMet", "StMet", 100, 0.0, 3000.0)


	return (plot_PT_b1, plot_ETA_b1, plot_PHI_b1, plot_PT_tau1, plot_ETA_tau1, plot_PHI_tau1, plot_PT_tau2, plot_ETA_tau2, plot_PHI_tau2, plot_PT_j1, plot_ETA_j1, plot_PHI_j1, plot_PT_j2, plot_ETA_j2, plot_PHI_j2, plot_DeltaR_tau1tau2, plot_DeltaR_j1j2, plot_DeltaPhi_tau1tau2, plot_DeltaPhi_j1j2, plot_DeltaPt_tau1tau2, plot_DeltaPt_j1j2, plot_M_Zp, plot_M_j1j2, plot_M_j1j2b1, plot_StMet)
