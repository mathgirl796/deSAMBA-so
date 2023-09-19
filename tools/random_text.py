from random import choice

with open('a.txt', 'w') as f:
    for i in range(100000):
        f.write(choice('acgtACGTacgtACGTacgtACGTacgtACGTacgtACGTacgtACGTacgtACGTacgtACGTacgtACGT\t\n\r !@#$%^&*()_+'))