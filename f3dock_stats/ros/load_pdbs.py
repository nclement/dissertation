import glob

bound_name='bound.out'
cmd.load('bound.out.pdb', bound_name)
cmd.select('cont_gold_r', '%s & n. ca & (chain R & resi 191+199+202-207+209-210+213+216-231+263-280+299-308+311+323-344+354-355+358-359+362+397-406+408-409+414-415)' % bound_name)
cmd.color('marine', '%s & chain R' % bound_name)
cmd.color('tv_green', '%s & chain L' % bound_name)
cmd.color('tv_red', 'cont_gold_r')

for i, fl in enumerate(glob.glob('1BKD.out.tog*.pdb')):
    cmd.load(fl, 'test_%d' % i)
    cmd.select('cont_test_r_%d' % i,'test_%d & n. ca & (chain R & resi 214+222+225-230+232-233+236+239-254+286-303+322-331+334+346-367+377-378+381-382+385+420-429+431-432+437-438)' % i)
    cmd.align('test_%d & cont_test_r_%d' % (i, i), '%s & cont_gold_r' % bound_name)
    cmd.color('marine', 'test_%d & chain R' % i)
    cmd.color('tv_green', 'test_%d & chain L' % i)
    cmd.color('tv_yellow', 'cont_test_r_%d' % i)

cmd.group('test_pdbs', 'test_*')
