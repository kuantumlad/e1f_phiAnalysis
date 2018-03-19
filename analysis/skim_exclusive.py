#!/usr/bin/env python 

import ROOT 
import sys

def calc_missing(particles):
    final_state = ROOT.TLorentzVector(0.0, 0.0, 5.498, 5.498+0.938)
    
    for p in particles:
        final_state -= p
        
    return final_state

def passes(event):
    particles = [event.electron, event.proton, event.kaon_pos]
    miss = calc_missing(particles).M()
    
    if miss > 0.45 and miss < 0.55:
        return True

    return False

class HistoPack(object):

    def __init__(self, name):
        self.name = name
        self._setup_histos()
    

    def _setup_histos(self):
        self.missing_mass_epx = ROOT.TH1F('missing_mass_epx' + str(self.name), 
                                          'missing_mass_epx' + str(self.name),
                                          200, 0.0, 2.0)

        self.missing_mass_epkx = ROOT.TH1F('missing_mass_epkx' + str(self.name), 
                                          'missing_mass_epkx' + str(self.name),
                                          200, 0.0, 1.0)

        self.im_kk_im_kp = ROOT.TH2F('im_kk_im_kp' + str(self.name), 
                                     'im_kk_im_kp' + str(self.name), 
                                     200, 0.7, 2.0,
                                     200, 1.4, 3.0)

    def fill(self, event):
        self.missing_mass_epx.Fill( calc_missing(list([event.electron, event.proton])).M() )
        self.missing_mass_epkx.Fill(  calc_missing(list([event.electron, event.proton, event.kaon_pos])).M() )

        im_kk = (event.kaon_pos + event.kaon_neg).M()
        im_kp = (event.kaon_pos + event.proton).M()
        self.im_kk_im_kp.Fill(im_kk, im_kp)

    def save(self, out_tfile):
        self.missing_mass_epx.Write()
        self.missing_mass_epkx.Write()
        self.im_kk_im_kp.Write()

if __name__ == "__main__":

    input_file = '../data/events.root'
    
    # setup reader
    reader = ROOT.TChain('events')
    reader.AddFile(input_file)

    # setup histo 
    histos = {}
    histos['all'] = HistoPack('all')
    histos['cut'] = HistoPack('cut')

    # tell how many events 
    print('Processing %s with %d entries.' % (input_file, reader.GetEntries()))
    
    # loop 
    for index, event in enumerate(reader):

        # the topology for ep->epk+X 
        if event.topology is 0:
            histos['all'].fill(event)
            
            if passes(event):
                histos['cut'].fill(event)

        # update the user 
        if index%10000 is 0:
            sys.stdout.write('\r Processed (%d/%d) events' % (index, reader.GetEntries()))
            sys.stdout.flush()


    out_file = ROOT.TFile.Open('out.root', 'recreate')

    histos['all'].save(out_file)
    histos['cut'].save(out_file)

    out_file.Close()

