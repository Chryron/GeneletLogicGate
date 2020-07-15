from biocrnpyler import *

class TranscriptionSwitch(Mechanism):
    """
    Reactions involved in this mechanism:
    
    Sw_OFF + A -> Sw_ON
    Sw_ON + I -> Sw_OFF + AI
    A + I -> AI
    Sw_OFF + RNAP <-> Sw_OFF.RNAP -> Sw_OFF + RNAP + P
    Sw_ON + RNAP <-> Sw_ON.RNAP -> Sw_ON + RNAP + P
    AI + RNAseH <-> AI.RNAseH -> A + RNAseH
    
    Optional reactions depending on inputs:
    
    Sw_OFF + A2 -> Sw_ON_2
    Sw_ON_2 + I2 -> Sw_OFF + AI_2
    A2 + I2 -> AI_2
    
    Allows upto 2 sets of activators and inhibitors per genelet
    """
    
    # Set the name and mechanism_type
    def __init__(self, name="transcription_switch", mechanism_type="transcription"):
        
            Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)
    
    def update_species(self, switch_off, transcript, activator, inhibitor, rnap, rnaseH, activator2 = None, inhibitor2 = None, **keywords):
        
        # Create appropriate complex species depending on whether or not 2nd set of activators and inhibitors are present
        
        
        if activator2 != None and inhibitor2 != None:
            A_I_complex2 = ComplexSpecies([inhibitor2, activator2],name = str(switch_off).replace("_OFF","")+"_AI_2")
            switch_on = ComplexSpecies([switch_off, activator], name = str(switch_off).replace("_OFF","_ON_1"))
            switch_on2 = ComplexSpecies([switch_off, activator2], name = str(switch_off).replace("_OFF","_ON_2"))
            self.switch_on2 = switch_on2
            self.A_I_complex2 = A_I_complex2
        else:     
            switch_on = ComplexSpecies([switch_off, activator], name = str(switch_off).replace("_OFF","_ON"))
        

        A_I_complex = ComplexSpecies([inhibitor, activator],name = str(switch_off).replace("_OFF","")+"_AI")
        self.switch_on = switch_on
        self.A_I_complex = A_I_complex
        
        # Return appropriate species depending on whether or not 2nd set of activators and inhibitors are present 
        
        if activator2 != None and inhibitor2 != None:
            return [switch_off, switch_on, switch_on2, transcript, activator, inhibitor, rnap, rnaseH, activator2, inhibitor2, A_I_complex, A_I_complex2,
                   ComplexSpecies([rnap, switch_on]), ComplexSpecies([rnap, switch_off]), ComplexSpecies([rnaseH, A_I_complex]),
                   ComplexSpecies([rnap, switch_on2]), ComplexSpecies([rnaseH, A_I_complex2])] 
        else:
            return [switch_off, switch_on, transcript, activator, inhibitor, rnap, rnaseH, A_I_complex,
                   ComplexSpecies([rnap, switch_on]), ComplexSpecies([rnap, switch_off]), ComplexSpecies([rnaseH, A_I_complex])] 
            
    
    def update_reactions(self, switch_off, transcript, activator, inhibitor, component, part_id, rnap, rnaseH, activator2 = None, inhibitor2 = None, **keywords):
        
        if activator2 != None and inhibitor2 != None:
            A_I_complex2 = self.A_I_complex2
            switch_on2 = self.switch_on2     
            
        switch_on = self.switch_on 
        A_I_complex = self.A_I_complex
        
        # Initialise reaction parameters
        
        ktx = component.get_parameter("ktx", part_id = part_id, mechanism = self)
        kleak = component.get_parameter("kleak", part_id = part_id, mechanism = self)
        kdeg = component.get_parameter("kdeg", part_id = part_id, mechanism = self)
        
        ku_tx = component.get_parameter("ku_tx", part_id = part_id, mechanism = self)
        ku_leak = component.get_parameter("ku_leak", part_id = part_id, mechanism = self)
        ku_deg = component.get_parameter("ku_deg", part_id = part_id, mechanism = self)
        
        kM_tx = component.get_parameter("kM_tx", part_id = part_id, mechanism = self)
        kM_leak = component.get_parameter("kM_leak", part_id = part_id, mechanism = self)
        kM_deg = component.get_parameter("kM_deg", part_id = part_id, mechanism = self)
        
        kon = component.get_parameter("kon", part_id = part_id, mechanism = self)
        koff = component.get_parameter("koff", part_id = part_id, mechanism = self)
        ka = component.get_parameter("ka", part_id = part_id, mechanism = self)
            
        kb_tx = (ku_tx + ktx) / kM_tx
        kb_leak = (ku_leak + kleak) / kM_leak
        kb_deg = (ku_deg + kdeg) / kM_deg
        
        #ku_tx = 0.1
        #ku_leak = 0.1
        #ku_deg = 0.1
        
        # Create reactions
        
        reaction_activation_1 = Reaction(inputs = [switch_off, activator], outputs = [switch_on], 
                            k = kon )
        
        reaction_deactivation_1 = Reaction(inputs = [switch_on, inhibitor], outputs = [switch_off, A_I_complex], 
                            k = koff )
        
        reaction_complex_1 = Reaction(inputs = [activator, inhibitor], outputs = [A_I_complex], k = ka)
        
        # Step 1 of transcription and A-I degradation reactions for first set of A/I
        
        reaction_tx_1_1 = Reaction(inputs = [switch_on, rnap], outputs = [ComplexSpecies([rnap, switch_on])], k = kb_tx, k_rev = ku_tx )
        
        reaction_leak_1 = Reaction(inputs = [switch_off, rnap], outputs = [ComplexSpecies([rnap, switch_off])], k = kb_leak, k_rev = ku_leak)
        
        reaction_deg_1_1 = Reaction(inputs = [A_I_complex, rnaseH], outputs = [ComplexSpecies([rnaseH, A_I_complex])], k = kb_deg, k_rev = ku_deg)
        
        # Step 2 of transcription and A-I degradation reactions for first set of A/I
        
        reaction_tx_2_1 = Reaction(inputs = [ComplexSpecies([rnap, switch_on])], outputs = [switch_on, transcript, rnap], k = ktx)
        
        reaction_leak_2 = Reaction(inputs = [ComplexSpecies([rnap, switch_off])], outputs = [switch_off, transcript, rnap], k = kleak)
        
        reaction_deg_2_1 = Reaction(inputs = [ComplexSpecies([rnaseH, A_I_complex])], outputs = [rnaseH, activator], k = kdeg)
        
        # Reactions if second set of activators and inhibitors are present
        
        if activator2 != None and inhibitor2 != None:
            
            reaction_activation_2 = Reaction(inputs = [switch_off, activator2], outputs = [switch_on2], 
                            k = kon )
        
            reaction_deactivation_2 = Reaction(inputs = [switch_on2, inhibitor2], outputs = [switch_off, A_I_complex2], 
                            k = koff )
        
            reaction_complex_2 = Reaction(inputs = [activator2, inhibitor2], outputs = [A_I_complex2], k = ka)
            
            # Step 1 of transcription and A-I degradation reactions for second set of A/I
    
            reaction_tx_1_2 = Reaction(inputs = [switch_on2, rnap], outputs = [ComplexSpecies([rnap, switch_on2])], k = kb_tx, k_rev = ku_tx )
            
            reaction_deg_1_2 = Reaction(inputs = [A_I_complex2, rnaseH], outputs = [ComplexSpecies([rnaseH, A_I_complex2])], k = kb_deg, k_rev = ku_deg)
            
            # Step 2 of transcription and A-I degradation reactions for second set of A/I
        
            reaction_tx_2_2 = Reaction(inputs = [ComplexSpecies([rnap, switch_on2])], outputs = [switch_on2, transcript, rnap], k = ktx)
        
            reaction_deg_2_2 = Reaction(inputs = [ComplexSpecies([rnaseH, A_I_complex2])], outputs = [rnaseH, activator2], k = kdeg)
            
            return [reaction_activation_1, reaction_deactivation_1, reaction_complex_1, reaction_tx_1_1, reaction_tx_2_1, reaction_leak_1, reaction_leak_2,
                reaction_deg_1_1, reaction_deg_2_1, reaction_activation_2, reaction_deactivation_2, reaction_complex_2, reaction_tx_1_2, reaction_tx_2_2,
                reaction_deg_1_2, reaction_deg_2_2]
            
        return [reaction_activation_1, reaction_deactivation_1, reaction_complex_1, reaction_tx_1_1, reaction_tx_2_1, reaction_leak_1, reaction_leak_2,
                reaction_deg_1_1, reaction_deg_2_1]
  


class Genelet(Promoter):
    """
    Genelet switch component using TranscriptionSwitch() mechanism
    Arguments: name, transcript, activator, inhibitor
    Optional Arguments: activator2, inhibitor2, rnap, rnaseH
    """
    def __init__(self, name, transcript, activator, inhibitor, rnap="RNAP", rnaseH="RNAseH", activator2 = None, inhibitor2 = None,  **keywords):
        
        # Set the Regulator
        # Component.set_species(species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species, string, or Component.
        
        self.activator = self.set_species(activator, material_type = "dna") 
        self.inhibitor = self.set_species(inhibitor, material_type = "rna")
        self.transcript = self.set_species(transcript, material_type = "rna")
        self.switch_off = self.set_species(str(name)+"_OFF")
        self.rnap = self.set_species(rnap, material_type = "protein")
        self.rnaseH = self.set_species(rnaseH, material_type = "protein")
        
        # Set second activator and inhibitors depending on whether they are present in the input
        
        if activator2 != None and inhibitor2 != None:
            
            self.activator2 = self.set_species(activator2, material_type = "dna") 
            self.inhibitor2 = self.set_species(inhibitor2, material_type = "rna")
        else:
            self.activator2 = None
            self.inhibitor2 = None
        
        custom_mechanisms = {"transcription": TranscriptionSwitch()} 
        
        Promoter.__init__(self, name = name, transcript = transcript, mechanisms = custom_mechanisms, **keywords)

    def update_species(self, **keywords):
        
        mech_tx = self.mechanisms["transcription"]
        
        species = [] 
        
        # Call update_species with correct arguments depending on whether the second set of activator and inhibitor are present
        
        if self.activator2 != None and self.inhibitor2 != None:
            species += mech_tx.update_species(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                          rnap = self.rnap, rnaseH = self.rnaseH, activator2 = self.activator2, inhibitor2 = self.inhibitor2)
        else:
            species += mech_tx.update_species(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                          rnap = self.rnap, rnaseH = self.rnaseH)
        
        return species

    def update_reactions(self, **keywords):
        mech_tx = self.mechanisms["transcription"]
        
        reactions = []
        
        # Call update_reactions with correct arguments depending on whether the second set of activator and inhibitor are present
        
        if self.activator2 != None and self.inhibitor2 != None:
            reactions += mech_tx.update_reactions(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                                  rnap = self.rnap, rnaseH = self.rnaseH, activator2 = self.activator2, inhibitor2 = self.inhibitor2,
                                                  component = self, part_id = "Genelet", **keywords)
        else:
            reactions += mech_tx.update_reactions(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                                  rnap = self.rnap, rnaseH = self.rnaseH, component = self, part_id = "Genelet", **keywords)
        return reactions
    
class Source(Promoter):
    """
    Genelet source component using Transcription_MM() mechanism
    Arguments: name, transcript
    Optional argument: rnap
    """
    def __init__(self, name, transcript, rnap="RNAP", **keywords):
        
        # Set the inouts
        # Component.set_species(species, material_type = None, attributes = None)
        # is a helper function that allows the input to be a Species, string, or Component
        
        self.dna = self.set_species(name)
        self.rnap = self.set_species(rnap, material_type = "protein")
        
        custom_mechanisms = {"transcription": Transcription_MM()}
        
        Promoter.__init__(self, name = name, transcript = transcript, mechanisms = custom_mechanisms, **keywords)

    def update_species(self, **keywords):
        
        mech_tx = self.mechanisms["transcription"]
        
        species = [] 
        species += mech_tx.update_species(dna = self.dna, transcript = self.transcript, rnap = self.rnap)
        
        return species

    def update_reactions(self, **keywords):
        mech_tx = self.mechanisms["transcription"]
        
        reactions = [] 
        reactions += mech_tx.update_reactions(dna = self.dna, transcript = self.transcript, 
                                              rnap = self.rnap, component = self, part_id = "Source", **keywords)
        return reactions    
    
    
def GeneletGate(name, out, on_1 = None, on_2 = None, off_1 = None, off_2 = None, typ = "AND"):
    """
    Function to initialise a unique modular Genelet Gate that can be added to a mixture. 
    Arguments: Name of gate (acts as a prefix for all species names involved in the gate)
               Type of gate ("AND" or "OR")
               Name of the ouput transcript of the gate
    Optional arguments: Activator and Inhibitor names for input genelet switches
    Output: List containing required Switch and Source Components
            Dictionary containing required concentrations of components
    """
    # User input error checking
    
    if type(name) != str:
        raise RuntimeError('AND gate name must be a string')
    if typ != "AND" and typ != "OR":
        raise RuntimeError('Gate type, typ, must be AND or OR')
        
    # Initialising transcriptional source concentration based on type of gate    
    
    if typ == "AND":
        #source = 164
        source = 170
    elif typ == "OR":
        #source = 130
        source = 102
    
    # Initialising components if optional arguments are not provided   
        
    if off_1 == None:
        off_1 = name + "_I1"
    if off_2 == None:
        off_2 = name + "_I2"
    if on_1 == None:
        on_1 = name + "_A1"
    if on_2 == None:
        on_2 = name + "_A2"
    
    # Logic gate component and initial condiction dictionary creation
    
    S1_off = Species(name + "_INP1")
    S2_off = Species(name + "_INP2")
    S3_off = Species(name + "_OUT")
    So1_on = Species(name + "_SOU")

    S1 = Genelet(S1_off, transcript = name + "_out_A", activator = on_1, inhibitor = off_1 )
    S2 = Genelet(S2_off, transcript = name + "_out_A", activator = on_2, inhibitor = off_2 )
    S3 = Genelet(S3_off, transcript = out, activator = name + "_out_A", inhibitor = name + "_out_I" )
    So1 = Source(So1_on, transcript = name + "_out_I")
    
    ic = {str(S1_off)+"_OFF": 500, "rna_"+on_1: 700, "rna_"+off_1: 200, str(S2_off)+"_OFF": 500, "rna_"+on_2: 700, "rna_"+off_2: 200, str(S3_off)+"_OFF": 500,
          "rna_"+name+"_out_A": 0, "rna_"+name+"_out_I": 0, str(So1_on):source, "protein_RNAP":100}
    
    return [S1,S2,S3,So1],ic
    