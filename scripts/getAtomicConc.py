import numpy as np

#Al:

AtomicConc = (0.92/26.981539) / ((0.92/26.981539) + (0.055/65.38) + (0.025/24.305))
print("Al atomic percent : %12.6f" % (AtomicConc))
AtomicConc = (0.055/65.38) / ((0.92/26.981539) + (0.055/65.38) + (0.025/24.305))
print("Zn atomic percent : %12.6f" % (AtomicConc))
AtomicConc = (0.025/24.305) / ((0.92/26.981539) + (0.055/65.38) + (0.025/24.305))
print("Mg atomic percent : %12.6f" % (AtomicConc))
