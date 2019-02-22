

massScale = 300.0
length=100
def temp(n):
    return massScale/(1+999/length * n)

file = open("sigmaTest.txt",'w')

for i in range(0,length+1):
    line = str(temp(i))
    for j in range(105):
        if j==49:
            line += " 10.0"
        elif j in [79,81,83,96]:
            line += " 0.5"
        elif j in [98,103]:
            line += " 0.0005"
        elif j in [42,45,47]:
            line += " 0.000000000001"
        else:
            line += " 0"
    if i!= length:
        line+="\n"
    file.write(line)

file.close()
