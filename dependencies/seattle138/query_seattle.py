from subprocess import Popen, PIPE
import sys

def download_seattle(classpath,input,output):
    cmd = ["java", "-classpath", classpath, "SubmitSeattleSeqAnnotationAutoJob138", input, output] 
    print("Started querying SeattleSeq138.. \n")
    with Popen(cmd, stdout=PIPE, bufsize=1, universal_newlines=True) as p:
        for line in p.stdout:
            print(line, end='') 
    print("Finished querying SeattleSeq138.. \n")
