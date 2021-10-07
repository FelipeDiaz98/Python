import ROOT
from ROOT import TLorentzVector, TH1F
import histos

def event_reconstruction(jets, bs, taus):

	global dijet_1
	global b_dijet_1
	global dijet_2
	global b_dijet_2
	global reconstructed_W1
	global reconstructed_W2
	#merged category of the first pair of hadronic jets coming from the top.
	global notMerged_1
	global partiallyMerged_1
	global fullyMerged_1
	#merged category of the second pair of hadronic jets coming from the top.
	global notMerged_2
	global partiallyMerged_2
	global fullyMerged_2
	#index of the first b-jet from the hadronic decaying top
	global index_b1
	#index of the second b-jet from the hadronic decaying top
	global index_b2
	#index of the two pairs of jets coming from the hadronic tops (pairs are always j1-j2 and j3-j4).
	global index_j1
	global index_j2
	global index_j3
	global index_j4
	#index of the two taus (pairs of jets tagged as coming from taus), does NOT include the energy from the neutrinos (MET).
	global index_tau1
	global index_tau2

	#global Delta_PT_b_alpha_b_beta

	best_Err = 999999.9
	MW = 80.379
	Mt = 172.76

	index_j_temp1 = 0
	index_j_temp2 = 0
	index_j_temp3 = 0
	index_j_temp4 = 0
	
	index_b_temp1 = 0
	index_b_temp2 = 0
	
	#loop over the first group of particles (j_1, j_2, b_1).
	for i in range(len(jets)):
		for j in range(len(jets)):
			if j > i:
				for k in range(len(bs)):
				
					#Loop over the second group of particles (j_3, j_4, b_2)
					for n in range(len(jets)):
						if ((n != i) and (n != j)):
							for m in range(len(jets)):
								if ((m > n) and (m != i) and (m != j)):
									for l in range(len(bs)):
										if (l != k):
							
											jtemp1 = jets[i]
											jtemp2 = jets[j]
											jtemp3 = jets[n]
											jtemp4 = jets[m]
											btemp1 = bs[k]
											btemp2 = bs[l]
												
											Err_1 = (abs((jtemp1 + jtemp2 + btemp1).M() - Mt))*100/Mt
											Err_2 = (abs((jtemp3 + jtemp4 + btemp2).M() - Mt))*100/Mt
											Err = Err_1 + Err_2
					
											#Selection criteria.
											if (Err < best_Err):
												best_Err = Err
												index_j_temp1 = i
												index_j_temp2 = j
												index_j_temp3 = n
												index_j_temp4 = m
												index_b_temp1 = k
												index_b_temp2 = l
					
	
	#Initialize the first group of particles.				
	index_j1 = index_j_temp1
	index_j2 = index_j_temp2	
	index_b1 = index_b_temp1
	
	dijet_1 = jets[index_j1] + jets[index_j2]		
	dr_dijet_1 = jets[index_j1].DeltaR(jets[index_j2])
	b_dijet_1 = jets[index_j1] + jets[index_j2] + bs[index_b1]
	dr_b_dijet_1 = dijet_1.DeltaR(bs[index_b1])
	
	#Initialize the second group of particles.
	index_j3 = index_j_temp3
	index_j4 = index_j_temp4	
	index_b2 = index_b_temp2
	
	dijet_2 = jets[index_j3] + jets[index_j4]		
	dr_dijet_2 = jets[index_j3].DeltaR(jets[index_j4])
	b_dijet_2 = jets[index_j3] + jets[index_j4] + bs[index_b2]
	dr_b_dijet_2 = dijet_2.DeltaR(bs[index_b2])
    		
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
		
	#Merged category for the second group of particles.	
	if (dr_dijet_2 > 0.8):

	 	notMerged_2 = True
	 	partiallyMerged_2 = False
	 	fullyMerged_2 = False


   	else:

		notMerged_2 = False
		
        	if (dr_b_dijet_2 > 1.0):
        
			partiallyMerged_2 = True
			fullyMerged_2 = False
			reconstructed_W2 = dijet_2

		else:

			partiallyMerged_2 = False
			fullyMerged_2 = True


	best_dPt = 99999999.9
	for i in range(len(taus)):
		for j in range(len(taus)):
			if (j > i):
			
				dPt = abs(taus[i].Pt() - taus[j].Pt())
				if (dPt < best_dPt):
				
					best_dPt = dPt
					index_tau1 = i
					index_tau2 = j
					
	return jets[index_j1], jets[index_j2], jets[index_j3], jets[index_j4], bs[index_b1], bs[index_b2], taus[index_tau1], taus[index_tau2]


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
jobs = [2,2,2,2]



#------------------ HISTOGRAMS ---------------------
c1 = ROOT.TCanvas("c1", "Titulo")    # ROOT canvas

#Creation of empty TH1F objects (empty ROOT histograms)
plots = histos.histos()   # plots is a list of TH1F objects


#------------- ITERATING THE FILES AND MAKING THE HISTOGRAMS ----------------

for n_signal, signal in enumerate(signals):

	f = ROOT.TFile(signal + ".root", "recreate")

	for ind in range(1, jobs[n_signal] + 1):

		directory = str("/disco4/SIMULACIONES/Angelica/" + signal + "/" + signal + "_" + str(ind) + "/Events/run_01/tag_1_delphes_events.root")
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

			# MET (neutrinos)
			Total_MET = 0
			EntryFromBranch_MET = File.MissingET.GetEntries()
			for j in range(EntryFromBranch_MET):
				MET = TLorentzVector()
				MET_PT, MET_Eta, MET_Phi, MET_M  = File.GetLeaf("MissingET.MET").GetValue(j), File.GetLeaf("MissingET.Eta").GetValue(j), File.GetLeaf("MissingET.Phi").GetValue(j), 0.0
				MET.SetPtEtaPhiM(MET_PT, MET_Eta, MET_Phi, MET_M)
				METs.append(MET)
				Total_MET += MET_PT
      
      			
      			#Checks if the event has the minimum theoretical expected number of particles (the count might be higher due to particle-detector effects).
			if (len(jets) >= 4 and len(bs) >= 2 and len(taus) >= 2):

				jets.sort(reverse = True, key=PT)     
				bs.sort(reverse = True, key=PT)
				taus.sort(reverse = True, key=PT)

				j1, j2, j3, j4, b1, b2, tau1, tau2 = event_reconstruction(jets, bs, taus)

				#Variables to be plotted.
				variables = [b1.Pt(), b1.Eta(), b1.Phi(), b2.Pt(), b2.Eta(), b2.Phi(), tau1.Pt(), tau1.Eta(), tau1.Phi(), tau2.Pt(), tau2.Eta(), tau2.Phi(), j1.Pt(), j1.Eta(), j1.Phi(), j2.Pt(), j2.Eta(), j2.Phi(), j3.Pt(), j3.Eta(), j3.Phi(), j4.Pt(), j4.Eta(), j4.Phi(), tau1.DeltaR(tau2), b1.DeltaR(b2), j1.DeltaR(j2), j1.DeltaR(j3), j1.DeltaR(j4), j2.DeltaR(j3), j2.DeltaR(j4), j3.DeltaR(j4), tau1.DeltaPhi(tau2), b1.DeltaPhi(b2), j1.DeltaPhi(j2), j1.DeltaPhi(j3), j1.DeltaPhi(j4), j2.DeltaPhi(j3), j2.DeltaPhi(j4), j3.DeltaPhi(j4), (tau1 + tau2).M() + Total_MET, (j1+j2).M(), (j1+j3).M(), (j1+j4).M(), (j2+j3).M(), (j2+j4).M(), (j3+j4).M(), (j1+j2+b1).M(), (j1+j3+b1).M(), (j1+j4+b1).M(), (j2+j3+b1).M(), (j2+j4+b1).M(), (j3+j4+b1).M(), (j1+j2+b2).M(), (j1+j3+b2).M(), (j1+j4+b2).M(), (j2+j3+b2).M(), (j2+j4+b2).M(), (j3+j4+b2).M(), j1.Pt() + j2.Pt() + j3.Pt() + j4.Pt() + b1.Pt() + b2.Pt() + tau1.Pt() + tau2.Pt() + Total_MET]

				for i in range(len(plots)):
			  	  histos_fill(plots[i], variables[i])


	  # Drawing in the histograms

	for plot in plots:
	  histos_Draw(plot)


	  # Updating the canvas
	c1.Update()


	  # Writing in the histograms
	for plot in plots:
	  histos_Write(plot)
	  

	  # Closing the ROOT file where the histos were saved.
	f.Close()

	  
	  # Reseting the TH1F objects for its use in the next signal or bkg file
	for plot in plots:
	  histos_Reset(plot)


