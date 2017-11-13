# create asci s4 class object
setClass('asci', 
         slots = list(
           scores = 'data.frame',
           Supp1_mmi = 'data.frame', 
           Supp1_OE = 'data.frame', 
           Supp2_OE = 'data.frame', 
           Supp3_OE = 'data.frame',
           taxa = 'character'
           )
)

