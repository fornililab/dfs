load compensatory_power.pdb
hide everything

create chain_A, chain A

hide ev
bg_color white
show cartoon, chain_A
center visi

# Spectrum interval defined by  minimum=2nd percentile and  maximum=98th percentile 
spectrum b, palette = blue_white_red, selection = chain_A, minimum = 23.3, maximum = 64.8
show surface, chain_A
# Selection of top 28% residues, P > 49.1
create predicted_compensatory_residues, (b > 49.1) and chain_A
hide everything, predicted_compensatory_residues
show spheres, predicted_compensatory_residues
color green, predicted_compensatory_residues
set transparency, 0.3
