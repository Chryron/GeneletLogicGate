from biocrnpyler import *

class TranscriptionSwitch(Mechanism):
    
    def __init__(self, name="transcription_switch", mechanism_type="transcription"):
        
            Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)
    
    def update_species(self, switch_off, transcript, activator, inhibitor, rnap, rnaseH, activator2 = None, inhibitor2 = None, **keywords):
        
        if activator2 != None and inhibitor2 != None:
            return [switch_off, transcript, activator, inhibitor, rnap, rnaseH, activator2, inhibitor2] 
        else:
            return [switch_off, transcript, activator, inhibitor, rnap, rnaseH] 
            
    
    def update_reactions(self, switch_off, transcript, activator, inhibitor, component, part_id, rnap, rnaseH, activator2 = None, inhibitor2 = None, **keywords):
        

        if activator2 != None and inhibitor2 != None:
            A_I_complex2 = ComplexSpecies([inhibitor2, activator2],name = str(switch_off).replace("_OFF","")+"_AI_2")
            switch_on = ComplexSpecies([switch_off, activator], name = str(switch_off).replace("_OFF","_ON_1"))
            switch_on2 = ComplexSpecies([switch_off, activator2], name = str(switch_off).replace("_OFF","_ON_2"))
        else:     
            switch_on = ComplexSpecies([switch_off, activator], name = str(switch_off).replace("_OFF","_ON"))
        A_I_complex = ComplexSpecies([inhibitor, activator],name = str(switch_off).replace("_OFF","")+"_AI")
        
        
        
        
        ktx = component.get_parameter("ktx", part_id = part_id, mechanism = self)
        kleak = component.get_parameter("kleak", part_id = part_id, mechanism = self)
        kdeg = component.get_parameter("kdeg", part_id = part_id, mechanism = self)
        
        kb_tx = component.get_parameter("kb_tx", part_id = part_id, mechanism = self)
        kb_leak = component.get_parameter("kb_leak", part_id = part_id, mechanism = self)
        kb_deg = component.get_parameter("kb_deg", part_id = part_id, mechanism = self)
        
        kM_tx = component.get_parameter("kM_tx", part_id = part_id, mechanism = self)
        kM_leak = component.get_parameter("kM_leak", part_id = part_id, mechanism = self)
        kM_deg = component.get_parameter("kM_deg", part_id = part_id, mechanism = self)
        
        kon = component.get_parameter("kon", part_id = part_id, mechanism = self)
        koff = component.get_parameter("koff", part_id = part_id, mechanism = self)
        ka = component.get_parameter("ka", part_id = part_id, mechanism = self)
            
        #ku_tx = (kb_tx + ktx) / kM_tx
        #ku_leak = (kb_leak + kleak) / kM_leak
        #ku_deg = (kb_deg + kdeg) / kM_deg
        
        ku_tx = 0.1
        ku_leak = 0.1
        ku_deg = 0.1
        
        reaction_activation_1 = Reaction(inputs = [switch_off, activator], outputs = [switch_on], 
                            k = kon )
        
        reaction_deactivation_1 = Reaction(inputs = [switch_on, inhibitor], outputs = [switch_off, A_I_complex], 
                            k = koff )
        
        reaction_complex_1 = Reaction(inputs = [activator, inhibitor], outputs = [A_I_complex], k = ka)
        
        reaction_tx_1_1 = Reaction(inputs = [switch_on, rnap], outputs = [ComplexSpecies([rnap, switch_on])], k = kb_tx, k_rev = ku_tx )
        
        reaction_leak_1 = Reaction(inputs = [switch_off, rnap], outputs = [ComplexSpecies([rnap, switch_off])], k = kb_leak, k_rev = ku_leak)
        
        reaction_deg_1_1 = Reaction(inputs = [A_I_complex, rnaseH], outputs = [ComplexSpecies([rnaseH, A_I_complex])], k = kb_deg, k_rev = ku_deg)
        
        reaction_tx_2_1 = Reaction(inputs = [ComplexSpecies([rnap, switch_on])], outputs = [switch_on, transcript, rnap], k = ktx)
        
        reaction_leak_2 = Reaction(inputs = [ComplexSpecies([rnap, switch_off])], outputs = [switch_off, transcript, rnap], k = kleak)
        
        reaction_deg_2_1 = Reaction(inputs = [ComplexSpecies([rnaseH, A_I_complex])], outputs = [rnaseH, activator], k = kdeg)
        
        
        if activator2 != None and inhibitor2 != None:
            
            reaction_activation_2 = Reaction(inputs = [switch_off, activator2], outputs = [switch_on2], 
                            k = kon )
        
            reaction_deactivation_2 = Reaction(inputs = [switch_on2, inhibitor2], outputs = [switch_off, A_I_complex2], 
                            k = koff )
        
            reaction_complex_2 = Reaction(inputs = [activator2, inhibitor2], outputs = [A_I_complex2], k = ka)
    
            reaction_tx_1_2 = Reaction(inputs = [switch_on2, rnap], outputs = [ComplexSpecies([rnap, switch_on2])], k = kb_tx, k_rev = ku_tx )
            
            reaction_deg_1_2 = Reaction(inputs = [A_I_complex2, rnaseH], outputs = [ComplexSpecies([rnaseH, A_I_complex2])], k = kb_deg, k_rev = ku_deg)
        
            reaction_tx_2_2 = Reaction(inputs = [ComplexSpecies([rnap, switch_on2])], outputs = [switch_on2, transcript, rnap], k = ktx)
        
            reaction_deg_2_2 = Reaction(inputs = [ComplexSpecies([rnaseH, A_I_complex2])], outputs = [rnaseH, activator2], k = kdeg)
            
            return [reaction_activation_1, reaction_deactivation_1, reaction_complex_1, reaction_tx_1_1, reaction_tx_2_1, reaction_leak_1, reaction_leak_2,
                reaction_deg_1_1, reaction_deg_2_1, reaction_activation_2, reaction_deactivation_2, reaction_complex_2, reaction_tx_1_2, reaction_tx_2_2,
                reaction_deg_1_2, reaction_deg_2_2]
            
        return [reaction_activation_1, reaction_deactivation_1, reaction_complex_1, reaction_tx_1_1, reaction_tx_2_1, reaction_leak_1, reaction_leak_2,
                reaction_deg_1_1, reaction_deg_2_1]
  


class Switch(Promoter):

    def __init__(self, name, transcript, activator, inhibitor, rnap="RNAP", rnaseH="RNAseH", activator2 = None, inhibitor2 = None,  **keywords):
        #Set the Regulator
        #Component.set_species(species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species, string, or Component.
        
        self.activator = self.set_species(activator, material_type = "dna") 
        self.inhibitor = self.set_species(inhibitor, material_type = "rna")
        self.transcript = self.set_species(transcript, material_type = "rna")
        self.switch_off = self.set_species(str(name)+"_OFF")
        self.rnap = self.set_species(rnap, material_type = "protein")
        self.rnaseH = self.set_species(rnaseH, material_type = "protein")
        
        if activator2 != None and inhibitor2 != None:
            
            self.activator2 = self.set_species(activator2, material_type = "dna") 
            self.inhibitor2 = self.set_species(inhibitor2, material_type = "rna")
        else:
            self.activator2 = None
            self.inhibitor2 = None
        
        #Mechanisms are inherited from the Mixture unless set specifically in self.default_mechanisms.
        custom_mechanisms = {"transcription": TranscriptionSwitch()}
        
        #Always call the superclass __init__() with **keywords passed through
        Promoter.__init__(self, name = name, transcript = transcript, mechanisms = custom_mechanisms, **keywords)

    def update_species(self, **keywords):
        #Mechanisms are stored in an automatically created dictionary: mechanism_type --> Mechanism Instance.
        mech_tx = self.mechanisms["transcription"]
        
        species = [] #A list of species must be returned
        if self.activator2 != None and self.inhibitor2 != None:
            species += mech_tx.update_species(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                          rnap = self.rnap, rnaseH = self.rnaseH, activator2 = self.activator2, inhibitor2 = self.inhibitor2)
        else:
            species += mech_tx.update_species(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                          rnap = self.rnap, rnaseH = self.rnaseH)
        
        return species

    def update_reactions(self, **keywords):
        mech_tx = self.mechanisms["transcription"]
        
        reactions = [] #a list of reactions must be returned
        if self.activator2 != None and self.inhibitor2 != None:
            reactions += mech_tx.update_reactions(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                                  rnap = self.rnap, rnaseH = self.rnaseH, activator2 = self.activator2, inhibitor2 = self.inhibitor2,
                                                  component = self, part_id = "Switch", **keywords)
        else:
            reactions += mech_tx.update_reactions(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                                  rnap = self.rnap, rnaseH = self.rnaseH, component = self, part_id = "Switch", **keywords)
        return reactions
    
class Source(Promoter):
    def __init__(self, name, transcript, rnap="RNAP", **keywords):
        #Set the Regulator
        #Component.set_species(species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species, string, or Component
        self.dna = self.set_species(name)
        self.rnap = self.set_species(rnap, material_type = "protein")
        #Mechanisms are inherited from the Mixture unless set specifically in self.default_mechanisms.
        custom_mechanisms = {"transcription": Transcription_MM()}
        
        #Always call the superclass __init__() with **keywords passed through
        Promoter.__init__(self, name = name, transcript = transcript, mechanisms = custom_mechanisms, **keywords)

    def update_species(self, **keywords):
        #Mechanisms are stored in an automatically created dictionary: mechanism_type --> Mechanism Instance.
        mech_tx = self.mechanisms["transcription"]
        
        species = [] #A list of species must be returned
        species += mech_tx.update_species(dna = self.dna, transcript = self.transcript, rnap = self.rnap)
        
        return species

    def update_reactions(self, **keywords):
        mech_tx = self.mechanisms["transcription"]
        
        reactions = [] #a list of reactions must be returned
        reactions += mech_tx.update_reactions(dna = self.dna, transcript = self.transcript, 
                                              rnap = self.rnap, component = self, part_id = "Source", **keywords)
        return reactions    
    
    
def GeneletGate(name, out, on_1 = None, on_2 = None, off_1 = None, off_2 = None, typ = "AND"):
    if type(name) != str:
        raise RuntimeError('AND gate name must be a string')
    if typ != "AND" and typ != "OR":
        raise RuntimeError('Gate type, typ, must be AND or OR')
        
    if typ == "AND":
        #source = 164
        source = 170
    elif typ == "OR":
        #source = 130
        source = 102
        
    if off_1 == None:
        off_1 = name + "_I1"
    if off_2 == None:
        off_2 = name + "_I2"
    if on_1 == None:
        on_1 = name + "_A1"
    if on_2 == None:
        on_2 = name + "_A2"
    
    
    
    S1_off = Species(name + "_INP1")
    S2_off = Species(name + "_INP2")
    S3_off = Species(name + "_OUT")
    So1_on = Species(name + "_SOU")

    S1 = Switch(S1_off, transcript = name + "_out_A", activator = on_1, inhibitor = off_1 )
    S2 = Switch(S2_off, transcript = name + "_out_A", activator = on_2, inhibitor = off_2 )
    S3 = Switch(S3_off, transcript = out, activator = name + "_out_A", inhibitor = name + "_out_I" )
    So1 = Source(So1_on, transcript = name + "_out_I")
    
    ic = {str(S1_off)+"_OFF": 500, "rna_"+on_1: 700, "rna_"+off_1: 200, str(S2_off)+"_OFF": 500, "rna_"+on_2: 700, "rna_"+off_2: 200, str(S3_off)+"_OFF": 500,
          "rna_"+name+"_out_A": 0, "rna_"+name+"_out_I": 0, str(So1_on):source, "protein_RNAP":100}
    
    return [S1,S2,S3,So1],ic