"""Code to work with motifs
"""


import numpy as np

class PWM(object):
    def __init__(self, weights, name=None, threshold=None):
        self.weights = weights
        self.name = name
        self.threshold = threshold

    @staticmethod
    def from_homer_motif(motif_file):
        with open(motif_file) as fp:
            header = fp.readline().strip().split('\t')
            name = header[1]
            threshold = float(header[2])
            weights = np.loadtxt(fp)

        return PWM(weights, name, threshold)

    @staticmethod
    def get_encode_pwms(motif_file):
        pwms = []

        with open(motif_file) as fp:
            line = fp.readline().strip()
            while True:
                if line == '':
                    break

                header = line.strip('>').strip()
                weights = []
                while True:
                    line = fp.readline()
                    if line == '' or line[0] == '>':
                        break
                    weights.append(map(float, line.split()))
                pwms.append(PWM(np.array(weights).transpose(1,0), header))

        return pwms

    @staticmethod

    def get_homer_pwms(motif_file):
        pwms = {}

        with open(motif_file) as fp:
            line = fp.readline().strip()
            while True:
                if line == '':
                    break

                header = line.strip('>').strip()
                motif_name = header.split()[1]
                weights = []
                while True:
                    line = fp.readline()
                    if line == '' or line[0] == '>':
                        break
                    weights.append(map(float, line.split()))

                pwms[motif_name] = PWM(np.array(weights).transpose(1,0), motif_name)

        return pwms

    @staticmethod
    def from_cisbp_motif(motif_file):
        name = os.path.basename(motif_file)
        with open(motif_file) as fp:
            _ = fp.readline()
            weights = np.loadtxt(fp)[:, 1:]
        return PWM(weights, name)

    def get_weights(self):
        return self.weights



def idx_to_string(val):
    """Convert index to base pair
    """
    idx_to_string = {0: "A",
                     1: "C",
                     2: "G",
                     3: "T"}

    return idx_to_string[val]

    

def get_consensus_from_weights(pwm_weights):
    """Given PWM weights get back a consensus sequence
    """
    sequence = []
    for idx in range(pwm_weights.shape[1]):
        sequence.append(idx_to_string(np.argmax(pwm_weights[:,idx])))
        
    return ''.join(sequence)
