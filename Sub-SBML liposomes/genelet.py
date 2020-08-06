from biocrnpyler import *

class TranscriptionSwitch(Mechanism):
    """
    Reactions involved in this mechanism:
    
    Sw_OFF + A -> Sw_ON
    Sw_ON + I -> Sw_OFF + AI
    A + I -> AI
    
    Optional reactions depending on inputs:
    
    Sw_OFF + A2 -> Sw_ON_2
    Sw_ON_2 + I2 -> Sw_OFF + AI_2
    A2 + I2 -> AI_2
    
    Allows upto 2 sets of activators and inhibitors per genelet
    """
    
    # Set the name and mechanism_type
    def __init__(self, name="transcription_switch", mechanism_type="transcription"):
        
            Mechanism.__init__(self, name=name, mechanism_type=mechanism_type)
    
    def update_species(self, switch_off, switch_on, transcript, activator, inhibitor, rnap, rnaseH, A_I_complex, switch_on2 = None, A_I_complex2 = None, activator2 = None,
                       inhibitor2 = None, **keywords):
        
        
        # Return appropriate species depending on whether or not 2nd set of activators and inhibitors are present 
        
        if activator2 != None and inhibitor2 != None:
            return [switch_off, switch_on, switch_on2, transcript, activator, inhibitor, rnap, rnaseH, activator2, inhibitor2, A_I_complex, A_I_complex2] 
        else:
            return [switch_off, switch_on, transcript, activator, inhibitor, rnap, rnaseH, A_I_complex] 
            
    
    def update_reactions(self, switch_off, switch_on, transcript, activator, inhibitor, rnap, rnaseH, A_I_complex, part_id, component = None, switch_on2 = None,
                         A_I_complex2 = None, activator2 = None, inhibitor2 = None, **keywords):
        
        # Initialise reaction parameters
        
        kon = component.get_parameter("kon", part_id = part_id, mechanism = self)
        koff = component.get_parameter("koff", part_id = part_id, mechanism = self)
        ka = component.get_parameter("ka", part_id = part_id, mechanism = self)
        
        # Create reactions
        
        reaction_activation_1 = Reaction(inputs = [switch_off, activator], outputs = [switch_on], 
                            k = kon )
        reaction_deactivation_1 = Reaction(inputs = [switch_on, inhibitor], outputs = [switch_off, A_I_complex], 
                            k = koff )
        reaction_complex_1 = Reaction(inputs = [activator, inhibitor], outputs = [A_I_complex], k = ka)
        
        # Reactions if second set of activators and inhibitors are present
        
        if activator2 != None and inhibitor2 != None:
            
            reaction_activation_2 = Reaction(inputs = [switch_off, activator2], outputs = [switch_on2], 
                            k = kon )
            reaction_deactivation_2 = Reaction(inputs = [switch_on2, inhibitor2], outputs = [switch_off, A_I_complex2], 
                            k = koff )
            reaction_complex_2 = Reaction(inputs = [activator2, inhibitor2], outputs = [A_I_complex2], k = ka)
            
            return [reaction_activation_1, reaction_deactivation_1, reaction_complex_1, reaction_activation_2, reaction_deactivation_2, reaction_complex_2]
            
        return [reaction_activation_1, reaction_deactivation_1, reaction_complex_1]
  


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
        
        
        if activator2 != None and inhibitor2 != None:
            self.activator2 = self.set_species(activator2, material_type = "dna") 
            self.inhibitor2 = self.set_species(inhibitor2, material_type = "rna")
            
            A_I_complex2 = ComplexSpecies([self.inhibitor2, self.activator2],name = str(name) + "_AI_2")
            switch_on = ComplexSpecies([self.switch_off, self.activator], name = str(name) + "_ON_1")
            switch_on2 = ComplexSpecies([self.switch_off, self.activator2], name = str(name) + "_ON_2")
            
            self.switch_on2 = switch_on2
            self.A_I_complex2 = A_I_complex2
        else:     
            switch_on = ComplexSpecies([self.switch_off, self.activator], name = str(name) + "_ON")
        
        A_I_complex = ComplexSpecies([self.inhibitor, self.activator], name = str(name) +"_AI")
        
        self.switch_on = switch_on
        self.A_I_complex = A_I_complex
        
        
        
        # Set second activator and inhibitors depending on whether they are present in the input
        
        if activator2 != None and inhibitor2 != None:
            
            self.activator2 = self.set_species(activator2, material_type = "dna") 
            self.inhibitor2 = self.set_species(inhibitor2, material_type = "rna")
        else:
            self.activator2 = None
            self.inhibitor2 = None
        
        custom_mechanisms = {"transcription": TranscriptionSwitch(), "catalysis": MichalisMentenCopy(), "degradation": MichalisMenten()} 
        
        Promoter.__init__(self, name = name, transcript = transcript, mechanisms = custom_mechanisms, **keywords)

    def update_species(self, **keywords):
        
        mech_tx = self.mechanisms["transcription"]
        mech_cat = self.mechanisms["catalysis"]
        mech_deg = self.mechanisms["degradation"]

        
        species = [] 
        
        # Call update_species with correct arguments depending on whether the second set of activator and inhibitor are present
        
        if self.activator2 != None and self.inhibitor2 != None:
            species += mech_tx.update_species(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                              rnap = self.rnap, rnaseH = self.rnaseH, activator2 = self.activator2, inhibitor2 = self.inhibitor2,
                                              switch_on = self.switch_on, A_I_complex = self.A_I_complex, switch_on2 = self.switch_on2, A_I_complex2 = self.A_I_complex2)
            
            species += mech_cat.update_species(Enzyme = self.rnap, Sub = self.switch_on2, Prod = self.transcript)
            species += mech_deg.update_species(Enzyme = self.rnaseH, Sub = self.A_I_complex2, Prod = self.activator2)
        
            
        else:
            species += mech_tx.update_species(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                          rnap = self.rnap, rnaseH = self.rnaseH, switch_on = self.switch_on, A_I_complex = self.A_I_complex)
            
        species += mech_cat.update_species(Enzyme = self.rnap, Sub = self.switch_on, Prod = self.transcript)
        species += mech_cat.update_species(Enzyme = self.rnap, Sub = self.switch_off, Prod = self.transcript)
        species += mech_deg.update_species(Enzyme = self.rnaseH, Sub = self.A_I_complex, Prod = self.activator)
        
        
        return species

    def update_reactions(self, **keywords):
        
        mech_tx = self.mechanisms["transcription"]
        mech_cat = self.mechanisms["catalysis"]
        mech_deg = self.mechanisms["degradation"]
        part_id = "Genelet"
        
        ktx = self.get_parameter("ktx", part_id = part_id)
        kleak = self.get_parameter("kleak", part_id = part_id)
        kdeg = self.get_parameter("kdeg", part_id = part_id)
        
        ku_tx = self.get_parameter("ku_tx", part_id = part_id)
        ku_leak = self.get_parameter("ku_leak", part_id = part_id)
        ku_deg = self.get_parameter("ku_deg", part_id = part_id)
        
        kM_tx = self.get_parameter("kM_tx", part_id = part_id)
        kM_leak = self.get_parameter("kM_leak", part_id = part_id)
        kM_deg = self.get_parameter("kM_deg", part_id = part_id)
            
        kb_tx = (ku_tx + ktx) / kM_tx
        kb_leak = (ku_leak + kleak) / kM_leak
        kb_deg = (ku_deg + kdeg) / kM_deg
        
        
        reactions = []
        
        # Call update_reactions with correct arguments depending on whether the second set of activator and inhibitor are present
        
        if self.activator2 != None and self.inhibitor2 != None:
            reactions += mech_tx.update_reactions(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                                  rnap = self.rnap, rnaseH = self.rnaseH, activator2 = self.activator2, inhibitor2 = self.inhibitor2,
                                                  switch_on = self.switch_on, A_I_complex = self.A_I_complex, switch_on2 = self.switch_on2, A_I_complex2 = self.A_I_complex2,
                                                  component = self, part_id = part_id, **keywords)
            reactions += mech_cat.update_reactions(Enzyme = self.rnap, Sub = self.switch_on2, Prod = self.transcript, kb = kb_tx, ku = ku_tx, kcat = ktx)
            reactions += mech_deg.update_reactions(Enzyme = self.rnaseH, Sub = self.A_I_complex2, Prod = self.activator2, kb = kb_deg, ku = ku_deg, kcat = kdeg)
        
        else:
            reactions += mech_tx.update_reactions(switch_off = self.switch_off, transcript = self.transcript, activator = self.activator, inhibitor = self.inhibitor, 
                                                  rnap = self.rnap, rnaseH = self.rnaseH, switch_on = self.switch_on, A_I_complex = self.A_I_complex,
                                                  component = self, part_id = part_id, **keywords)
            
        reactions += mech_cat.update_reactions(Enzyme = self.rnap, Sub = self.switch_on, Prod = self.transcript, kb = kb_tx, ku = ku_tx, kcat = ktx)
        reactions += mech_cat.update_reactions(Enzyme = self.rnap, Sub = self.switch_off, Prod = self.transcript, kb = kb_leak, ku = ku_leak, kcat = kleak)
        reactions += mech_deg.update_reactions(Enzyme = self.rnaseH, Sub = self.A_I_complex, Prod = self.activator, kb = kb_deg, ku = ku_deg, kcat = kdeg)
        
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
    