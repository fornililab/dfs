load compensatory_power.pdb
remove chain B

select chain_A, chain A
select dna, chain E or chain F
disable chain_A
disable dna

hide ev
bg_color white
show cartoon, chain_A
show cartoon, dna

spectrum b, palette=blue_red, selection=chain_A
show surface, chain_A
set transparency, 0.3