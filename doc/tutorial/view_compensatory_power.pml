load compensatory_power.pdb
remove chain B

select chain_A, chain A
disable chain_A

hide ev
bg_color white
show cartoon, chain_A
center visi

spectrum b, palette=blue_red, selection=chain_A
show surface, chain_A
select predicted_compensatory_residues, b>4.90
color green, predicted_compensatory_residues
set transparency, 0.3
