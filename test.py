text = """
                DG   .D    ;:           
              E   DKEEE EE.   E         
               GEEEE.;EE E     D        
           EEEEEEEEf  E.EEEEEEE         
           EEEEK.  EK  K   K   K        
           EE            ;EEEED         
         EEE                EEE         
          E                    E        
         E                      E       
       ED                       KE      
      E E                        E      
      EE                         E      
     E                           E      
     EE                          EE     
    D K                         KEE     
    EDt   EDEL                   E.     
    EE   .E   E                  EEE  E 
  E EE   E    E                  EED   D
 E  EE   E         E             iEEE  E
 j  EE   E         EE      E K.   EEE  E
E   EK   D         EE     E   E   EEE  D
E  EEE                        E   KEE E 
 E E.E                            fEEE  
  EEEE                            ;KE   
     E               E            EEE   
     j         .     E            EEE   
      K        E     E            . E   
      E        i    E            E  E   
                E  .D           E       
       E         EE            .D       
        E                     .D        
        E                     E         
         E                   E          
          E                 E           
           E               E            
             K           E              
               EE      tE               
                  .KEEE                 
"""
padded = ''
for c in text:
    padded += ' ' + c

with open('sensei.txt', 'w') as f:
    f.write(padded)