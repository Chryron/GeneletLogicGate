from sympy import Symbol,sympify
from libsbml import *
import sys
import numpy as np

def sbml_to_ode2(filename):
    '''A function that takes in an SBML file and returns x,f,P,params_values.
    x is a list of species written as Sympy objects
    f is a list of functions written as Sympy objects
    P is a list of parameters written as Sympy objects
    params_values is a list of parameter values, in the same order as P
    x_init is a list of initial conditions, in the same order as x'''

    # Get the sbml file, check for errors, and perform conversions
    doc = readSBMLFromFile(filename)
    if doc.getNumErrors(LIBSBML_SEV_FATAL):
        print('Encountered serious errors while reading file')
        print(doc.getErrorLog().toString())
        sys.exit(1)
    
    doc.getErrorLog().clearLog()
    
    # Convert local params to global params
    props = ConversionProperties()
    props.addOption("promoteLocalParameters", True)
  
    if doc.convert(props) != LIBSBML_OPERATION_SUCCESS: 
        print('The document could not be converted')
        print(doc.getErrorLog().toString())
    
    # Expand initial assignments
    props = ConversionProperties()
    props.addOption("expandInitialAssignments", True)
  
    if doc.convert(props) != LIBSBML_OPERATION_SUCCESS: 
        print('The document could not be converted')
        print(doc.getErrorLog().toString())
    
    # Expand functions definitions
    props = ConversionProperties()
    props.addOption("expandFunctionDefinitions", True)
  
    if doc.convert(props) != LIBSBML_OPERATION_SUCCESS: 
        print('The document could not be converted')
        print(doc.getErrorLog().toString())
    
    # Get model and define importnat lists, dictionaries
    mod = doc.getModel()
    x = []
    x_init = []
    P = []
    params_values = []
    reactions = {}
 
    # Append species symbol to 'x' and append initial amount/concentration to x_init
    # x[i] corresponds to x_init[i]
    for i in range(mod.getNumSpecies()): 
        species = mod.getSpecies(i)
        x.append(Symbol(species.getId()))
        if species.isSetInitialConcentration():
            x_init.append(species.getInitialConcentration())
        elif species.isSetInitialAmount():
            x_init.append(species.getInitialAmount())
        else:
            x_init.append(0)

    # Append parameter symbol to 'P' and parameter values to 'params_values'
    for i in range(mod.getNumParameters()): 
        params = mod.getParameter(i)
        params_values.append(params.getValue())
        P.append(Symbol(params.getId()))

    
    
    
    # Get kinetic formula for each reaction, store in dictionary 'reactions'
    for i in range(mod.getNumReactions()): 
        reaction = mod.getReaction(i)
        kinetics = reaction.getKineticLaw()
        reactions[reaction.getId()] = sympify(kinetics.getFormula())

    # Define f
    f = [0] * len(x)
    
    # Loop to define functions in 'f'
    for i in range(mod.getNumReactions()): 
        reaction = mod.getReaction(i)
        # subtract reactant kinetic formula
        for j in range(reaction.getNumReactants()):
            ref = reaction.getReactant(j)
            species = sympify(mod.getSpecies(ref.getSpecies()).getId())
            curr_index = x.index(species)
            # Check stoichiometry
            if ref.getStoichiometry() == 1.0:  
                f[curr_index] += -reactions[reaction.getId()]
            else:
                f[curr_index] += -reactions[reaction.getId()] * ref.getStoichiometry()
        # add product kinetic formula
        for j in range(reaction.getNumProducts()):
            ref = reaction.getProduct(j)
            species = sympify(mod.getSpecies(ref.getSpecies()).getId())
            curr_index = x.index(species)
            # Check stoichiometry
            if ref.getStoichiometry() == 1.0: 
                f[curr_index] += +reactions[reaction.getId()]
            else:
                f[curr_index] += +reactions[reaction.getId()] * ref.getStoichiometry()


    return x,f,P,params_values,x_init