import reaktoro as rkt
import numpy as np

# add constraint for total amount of aqueous phase
def addSConstraint(specs): 
    idx_S = specs.addInput("S")

    S_constraint = rkt.ConstraintEquation()
    S_constraint.id = "S"
    
    def SFunction(props, w): 
        aq_props = rkt.AqueousProps(props)
        model_S_conc = aq_props.elementMolality("S") 
        return w[idx_S] - model_S_conc
        
    S_constraint.fn = SFunction   
    specs.addConstraint(S_constraint)


# create object to store constraint info
class Constraint: 
    def __init__(self, name, value, units): 
        self.name = name
        self.value = value
        self.units = units


# equilibrate state given system and series of constraints (i.e., conditions)
class ConstrainedEquilibration: 

    def __init__(self, system): 
        self.system = system
        self.specs = rkt.EquilibriumSpecs(system)
        self.constraints = []
        self.initial_state = rkt.ChemicalState(system)

    def addConstraint(self, name, value, units=""): 
        constraint = Constraint(name, value, units)
        self.constraints.append(constraint)

        if name == "temperature": 
            self.specs.temperature()   
        elif name == "pressure": 
            self.specs.pressure()
        elif name == "fO2": 
            self.specs.fugacity("O2(g)")
        elif name == "pH": 
            self.specs.pH()
        elif name == "S": 
            addSConstraint(self.specs)
            self.specs.openTo("S") # explicit titrant
        
        return value
    
    def initialize(self, substance, amount, units): 
        self.initial_state.add(substance, amount, units)

    def equilibrate(self): 
        self.solver = rkt.EquilibriumSolver(self.specs)
        self.conditions = rkt.EquilibriumConditions(self.specs)

        constraint_values = []
        for constraint in self.constraints: 
            constraint_values.append(constraint.value)
        
        mesh_values = np.meshgrid(*constraint_values)

        flattened_values = []
        for dim in mesh_values: 
            flattened_values.append(dim.flatten())

        self.states = []
        self.equilibrium_status = []

        for condition in zip(*flattened_values): 
            state = rkt.ChemicalState(self.initial_state)

            for value, constraint in zip(condition, self.constraints): 
                if constraint.name == "temperature": 
                    self.conditions.temperature(value, constraint.units)
                elif constraint.name == "pressure": 
                    self.conditions.pressure(value, constraint.units)   
                elif constraint.name == "fO2": 
                    self.conditions.fugacity("O2(g)", value, constraint.units)
                elif constraint.name == "pH": 
                    self.conditions.pH(value)          
                elif constraint.name == "S": 
                    self.conditions.set(constraint.name, value) # [molal]

            state = rkt.ChemicalState(state)   
            
            result = self.solver.solve(state, self.conditions)

            self.states.append(state)
            self.equilibrium_status.append("Successful computation!" if result.optima.succeeded else "Computation has failed!")

        
        return self.states, self.equilibrium_status