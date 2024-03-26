# Create model using single template

from modeller import *
from modeller.automodel import *
import os

log.verbose()
env = Environ()

# Read in HETATM records from template PDBs
env.io.hetatm = True

code = '3m01'
seg = 'A'
code_tmpl = '2ong'+seg
# Note that the gaps are numbered as in model (so if you have initial offset, that
# must not be included in the numbering)
# Additional gaps added here based on template protein
gap = str()
gaps = str()
gaps = [(4, 13), (523, 529)]#, (566, 566)] # Original gap numbering 

selection_str = str()
# Single residue gaps
for i in range(len(gap)):
    selection_str += "self.residues[" + "'" + str(gap[i]) + ":" + seg + "'" + "]"
    if i < len(gap)-1:
        selection_str += ","

# Multiple residues gaps
for j in range(len(gaps)):
    if len(selection_str) > 0 and j == 0:
        selection_str += ","
    s = "'" + str(gaps[j][0]-3) + ":" + seg + "'" + ","
    selection_str += "self.residue_range(" + s
    s = "'" + str(gaps[j][1]+1) + ":" + seg + "'"
    selection_str += s + ")"
    if j < len(gaps)-1:
        selection_str += ","

# Define new class to define range that is flexible
class MyModel(AutoModel):
    def select_atoms(self):
        return eval('Selection(' + selection_str + ')')

a = MyModel(env, alnfile=code+'-'+code_tmpl+'-ligand.ali',
              knowns=code_tmpl, sequence=code,
              inifile=code+'_loop_model.pdb',
              assess_methods=(assess.DOPE,
                              #soap_protein_od.Scorer(),
                              assess.GA341))
a.starting_model = 1
a.ending_model = 50   # Should be bigger number (e.g., 50)
a.md_level = refine.very_slow  # very_fast, fast, slow, very_slow, slow_large
a.max_var_iterations = 500     # Number of CG steps

a.make()

# Get a list of all successfully built models from a.outputs
ok_models = [x for x in a.outputs if x['failure'] is None]

# Rank the models by DOPE score
key = 'DOPE score'
ok_models.sort(key=lambda a: a[key])

# Get top model
m = ok_models[0]
print("Top model: %s (DOPE score %.3f)" % (m['name'], m[key]))
command = 'cp ' + m['name'] + ' ' + code + '_final_model_ligand.pdb'
os.system(command)  
print(command)
#m.write(file=code+'_final_model.pdb')   # Write lowest energy model

