import ROOT
from ROOT import TLorentzVector, TH1F
import histos


def jet_reconstruction(jets, bs):

	global dijet_1
	global b_dijet_1
	global reconstructed_W_1
	#merged category of the first pair of hadronic jets coming from the top.
	global notMerged_1
	global partiallyMerged_1
	global fullyMerged_1
	#index of the first b-jet from the hadronic decaying top
	global index_b1
	#index of the pair of jets coming from the hadronic tops (pairs are always j1-j2 and j3-j4).
	global index_j1
	global index_j2

	#global Delta_PT_b_alpha_b_beta

	best_Err = 999999.9
	MW = 80.379
	Mt = 172.76

	index_j_temp1 = 0
	index_j_temp2 = 0
	
	index_b_temp1 = 0
	
	#loop over the first group of particles (j_1, j_2, b_1).
	for i in range(len(jets)):
		for j in range(i + 1, len(jets)):
				for k in range(len(bs)):
					
					jtemp1 = jets[i]
					jtemp2 = jets[j]
					btemp1 = bs[k]
						
					Err = (abs((jtemp1 + jtemp2 + btemp1).M() - Mt))*100/Mt

					#Selection criteria.
					if (Err < best_Err):
						best_Err = Err
						index_j_temp1 = i
						index_j_temp2 = j
						index_b_temp1 = k
					
	
	#Initialize the first group of particles.				
	index_j1 = index_j_temp1
	index_j2 = index_j_temp2	
	index_b1 = index_b_temp1
	
	dijet_1 = jets[index_j1] + jets[index_j2]		
	dr_dijet_1 = jets[index_j1].DeltaR(jets[index_j2])
	b_dijet_1 = jets[index_j1] + jets[index_j2] + bs[index_b1]
	dr_b_dijet_1 = dijet_1.DeltaR(bs[index_b1])
	
    		
	#Merged category for the first group of particles.
	if (dr_dijet_1 > 0.8):

	 	notMerged_1 = True
	 	partiallyMerged_1 = False
	 	fullyMerged_1 = False


   	else:

		notMerged_1 = False
		
        	if (dr_b_dijet_1 > 1.0):
        
			partiallyMerged_1 = True
			fullyMerged_1 = False
			reconstructed_W1 = dijet_1

		else:

			partiallyMerged_1 = False
			fullyMerged_1 = True
					
	return jets[index_j1], jets[index_j2], bs[index_b1]


def tau_reconstruction(taus):

	#index of the two taus (pairs of jets tagged as coming from taus), does NOT include the energy from the neutrinos (MET).
	global index_tau1
	global index_tau2
	
	best_dPt = 99999999.9
	for i in range(len(taus)):
		for j in range(i + 1, len(taus)):
			
				dPt = abs(taus[i].Pt() - taus[j].Pt())
				if (dPt < best_dPt):
				
					best_dPt = dPt
					index_tau1 = i
					index_tau2 = j
					
	return taus[index_tau1], taus[index_tau2]
	

def cross_cleaning(jets, bs, taus):

	particles = jets + bs
	temp = []
	for i in range(len(taus)):
		
		crossed = False
		for j in range(len(particles)):

			if (taus[i].DeltaR(particles[j]) < 0.3):
				crossed = True

		if (crossed == False):
			
			temp.append(taus[i])

	return temp


def PT(TLV):
    return TLV.Pt()
    
def histos_fill(plot, variable):
  return plot.Fill(variable)


def histos_Draw(plot):
    return plot.Draw('HISTOS')  


def histos_Write(plot):
    return plot.Write()


def histos_Reset(plot):
    return plot.Reset('ICESM') 


signals = ["ttbarh", "ttbarZ", "Zprime_tata_350", "Zprime_tata_1500"]
jobs = [20,20,20,20]



#------------------ HISTOGRAMS ---------------------
c1 = ROOT.TCanvas("c1", "Titulo")    # ROOT canvas

#Creation of empty TH1F objects (empty ROOT histograms)
plot_tausmalos = TH1F("tausmalos", "tausmalos", 2, 0.0, 3.0)
plot_P_tau2 = TH1F("P_tau2", "P_tau2", 100, 0.0, 1000.0)

plots = histos.histos()   # plots is a list of TH1F objects

plot_PT_tausmenores = TH1F("PT_tausmenores", "PT_tausmenores", 100, 0.0, 1000.0)
plot_PT_tausmayores = TH1F("PT_tausmayores", "PT_tausmayores", 100, 0.0, 1000.0)
plot_ETA_tausmenores = TH1F("ETA_tausmenores", "ETA_tausmenores", 100, -5, 5)
plot_ETA_tausmayores = TH1F("ETA_tausmayores", "ETA_tausmayores", 100, -5, 5)
plot_PHI_tausmenores = TH1F("PHI_tausmenores", "PHI_tausmenores", 100, -4, 4)
plot_PHI_tausmayores = TH1F("PHI_tausmayores", "PHI_tausmayores", 100, -4, 4)

plot_DeltaR_tausmenores = TH1F("DeltaR_tausmenores", "DeltaR_tausmenores", 100, -4, 4)
plot_DeltaR_tausmayores = TH1F("DeltaR_tausmayores", "DeltaR_tausmayores", 100, -4, 4)
plot_DeltaPhi_tausmenores = TH1F("DeltaPhi_tausmenores", "DeltaPhi_tausmenores", 100, -4, 4)
plot_DeltaPhi_tausmayores = TH1F("DeltaPhi_tausmayores", "DeltaPhi_tausmayores", 100, -4, 4)


#------------- ITERATING THE FILES AND MAKING THE HISTOGRAMS ----------------

for n_signal, signal in enumerate(signals):

	f = ROOT.TFile(signal + ".root", "recreate")

	for ind in range(1, jobs[n_signal] + 1):

		directory = str("/disco4/SIMULACIONES/Liliana/" + signal + "/" + signal + "_" + str(ind) + "/Events/run_01/tag_1_delphes_events.root")
		File = ROOT.TChain("Delphes;1")
		File.Add(directory)
		Number = File.GetEntries()

		print("Signal: " + signal + "_" + str(ind))

		for i in range(Number):
			Entry = File.GetEntry(i)

			#Initializes particles lists.
			jets = []
			bs = []
			METs = []
			taus = []
			
			tausbad = []

			EntryFromBranch_j = File.Jet.GetEntries()
			for j in range(EntryFromBranch_j):

				BTag = File.GetLeaf("Jet.BTag").GetValue(j)
				TauTag = File.GetLeaf("Jet.TauTag").GetValue(j)

				#searches for jets.
				if (BTag != 1 and TauTag != 1):
					jet = TLorentzVector()
					jet_PT, jet_Eta, jet_Phi, jet_M  = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
					jet.SetPtEtaPhiM(jet_PT, jet_Eta, jet_Phi, jet_M)
					jets.append(jet)

				#searches for b_jets.
				elif (BTag == 1 and TauTag != 1):
					bjet = TLorentzVector()
					bjet_PT, bjet_Eta, bjet_Phi, bjet_M  = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j)
					bjet.SetPtEtaPhiM(bjet_PT, bjet_Eta, bjet_Phi, bjet_M)
					bs.append(bjet)

				#searches for taus.
				elif (TauTag == 1 and BTag != 1):
					tau = TLorentzVector()
					tau_PT, tau_Eta, tau_Phi, tau_M, tau_Charge = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j),  File.GetLeaf("Jet.Charge").GetValue(j)
					tau.SetPtEtaPhiM(tau_PT, tau_Eta, tau_Phi, tau_M)
					taus.append(tau)
					
				#searches for jets with both tags, bad taus (tauTag and bTag).
				elif (TauTag == 1 and BTag == 1):
					tau = TLorentzVector()
					tau_PT, tau_Eta, tau_Phi, tau_M, tau_Charge = File.GetLeaf("Jet.PT").GetValue(j), File.GetLeaf("Jet.Eta").GetValue(j), File.GetLeaf("Jet.Phi").GetValue(j), File.GetLeaf("Jet.Mass").GetValue(j),  File.GetLeaf("Jet.Charge").GetValue(j)
					tau.SetPtEtaPhiM(tau_PT, tau_Eta, tau_Phi, tau_M)
					tausbad.append(tau)

			# MET (neutrinos)
			Total_MET = 0
			EntryFromBranch_MET = File.MissingET.GetEntries()
			for j in range(EntryFromBranch_MET):
				MET = TLorentzVector()
				MET_PT, MET_Eta, MET_Phi, MET_M  = File.GetLeaf("MissingET.MET").GetValue(j), File.GetLeaf("MissingET.Eta").GetValue(j), File.GetLeaf("MissingET.Phi").GetValue(j), 0.0
				MET.SetPtEtaPhiM(MET_PT, MET_Eta, MET_Phi, MET_M)
				METs.append(MET)
				Total_MET += MET_PT
				
				
			#print(len(taus), len(tausbad))	

			plot_tausmalos.AddBinContent(1, len(taus))
			
			plot_tausmalos.AddBinContent(2, len(tausbad))

			a = plot_tausmalos.GetXaxis()
			a.SetBinLabel(1,"Real taus")
			a.SetBinLabel(2,"Fake taus")
      
      			
      			#Checks if the event has the minimum theoretical expected number of particles (the count might be higher due to particle-detector effects).
			if (len(jets) >= 2 and len(bs) >= 2 and len(taus) >= 2):

				jets.sort(reverse = True, key=PT)     
				bs.sort(reverse = True, key=PT)
				taus.sort(reverse = True, key=PT)

				j1, j2, b1 = jet_reconstruction(jets, bs)
				realjets = [j1, j2]
				realbs = [b1]
				
				#Makes the cross cleaning, only saving taus with dR not less than 0.3 to any other jet.
				realtaus = cross_cleaning(realjets, realbs, taus)

				#Checks if there are at least two real taus in the event (taus that passed the cross cleaning).
				if (len(realtaus) >= 2):

					tau1, tau2 = tau_reconstruction(realtaus)
					
					#Pt requirement of the different jets (20 for taus is LHC minimum).
					if ((j1.Pt() > 30) and (j2.Pt() > 30) and (b1.Pt() > 30) and (tau1.Pt() > 20) and (tau2.Pt() > 20)):

						#LHC required eta for taus.				
						if ((abs(tau1.Eta()) < 2.3) and (abs(tau2.Eta()) < 2.3)):

							#Minimum possible separation between both taus.
							if (tau1.DeltaR(tau2) > 0.3):

								#Variables to be plotted.
								variables = [b1.Pt(), b1.Eta(), b1.Phi(), b2.Pt(), b2.Eta(), b2.Phi(), tau1.Pt(), tau1.Eta(), tau1.Phi(), tau2.Pt(), tau2.Eta(), tau2.Phi(), j1.Pt(), j1.Eta(), j1.Phi(), j2.Pt(), j2.Eta(), j2.Phi(), tau1.DeltaR(tau2), b1.DeltaR(b2), j1.DeltaR(j2), tau1.DeltaPhi(tau2), b1.DeltaPhi(b2), j1.DeltaPhi(j2), abs(tau1.Pt() - tau2.Pt()), abs(b1.Pt() - b2.Pt()), abs(j1.Pt() - j2.Pt()), (tau1 + tau2).M() + Total_MET, (j1+j2).M(), (j1+j2+b1).M(), (j1+j2+b2).M(), j1.Pt() + j2.Pt() + b1.Pt() + b2.Pt() + tau1.Pt() + tau2.Pt() + Total_MET]

								for i in range(len(plots)):
							  		histos_fill(plots[i], variables[i])
							  		
							  	histos_fill(plot_P_tau2, tau2.P())

								if (tau2.Pt() < 80.0):
									histos_fill(plot_PT_tausmenores, tau2.Pt())
									histos_fill(plot_ETA_tausmenores, tau2.Eta())
									histos_fill(plot_PHI_tausmenores, tau2.Phi())

									histos_fill(plot_DeltaR_tausmenores, tau1.DeltaR(tau2))
									histos_fill(plot_DeltaPhi_tausmenores, tau1.DeltaPhi(tau2))
								
								if (tau2.Pt() > 80.0):
									histos_fill(plot_PT_tausmayores, tau2.Pt())
									histos_fill(plot_ETA_tausmayores, tau2.Eta())
									histos_fill(plot_PHI_tausmayores, tau2.Phi())

									histos_fill(plot_DeltaR_tausmayores, tau1.DeltaR(tau2))
									histos_fill(plot_DeltaPhi_tausmayores, tau1.DeltaPhi(tau2))


	# Drawing in the histograms
	plot_tausmalos.Draw('HISTOS')
	plot_P_tau2.Draw('HISTOS')

	for plot in plots:
		histos_Draw(plot)
		
	plot_PT_tausmenores.Draw('HISTOS')
	plot_PT_tausmayores.Draw('HISTOS')
	plot_ETA_tausmenores.Draw('HISTOS')
	plot_ETA_tausmayores.Draw('HISTOS')
	plot_PHI_tausmenores.Draw('HISTOS')
	plot_PHI_tausmayores.Draw('HISTOS')

	plot_DeltaR_tausmenores.Draw('HISTOS')
	plot_DeltaR_tausmayores.Draw('HISTOS')
	plot_DeltaPhi_tausmenores.Draw('HISTOS')
	plot_DeltaPhi_tausmayores.Draw('HISTOS')


	# Updating the canvas
	c1.Update()


	# Writing in the histograms
	plot_tausmalos.Write()
	plot_P_tau2.Write()
	
	for plot in plots:
		histos_Write(plot)
	  
	plot_PT_tausmenores.Write()
	plot_PT_tausmayores.Write()
	plot_ETA_tausmenores.Write()
	plot_ETA_tausmayores.Write()
	plot_PHI_tausmenores.Write()
	plot_PHI_tausmayores.Write()

	plot_DeltaR_tausmenores.Write()
	plot_DeltaR_tausmayores.Write()
	plot_DeltaPhi_tausmenores.Write()
	plot_DeltaPhi_tausmayores.Write()
	  

	# Closing the ROOT file where the histos were saved.
	f.Close()

	  
	# Reseting the TH1F objects for its use in the next signal or bkg file
	plot_tausmalos.Reset('ICESM')
	plot_P_tau2.Reset('ICESM')
	
	for plot in plots:
		histos_Reset(plot)
		
	plot_PT_tausmenores.Reset('ICESM')
	plot_PT_tausmayores.Reset('ICESM')
	plot_ETA_tausmenores.Reset('ICESM')
	plot_ETA_tausmayores.Reset('ICESM')
	plot_PHI_tausmenores.Reset('ICESM')
	plot_PHI_tausmayores.Reset('ICESM')

	plot_DeltaR_tausmenores.Reset('ICESM')
	plot_DeltaR_tausmayores.Reset('ICESM')
	plot_DeltaPhi_tausmenores.Reset('ICESM')
	plot_DeltaPhi_tausmayores.Reset('ICESM')


