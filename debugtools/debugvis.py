# debugvis.py
import io
import base64
import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import lldb
import debugger

def show():

    image_bytes = io.BytesIO()
    plt.savefig(image_bytes, format='png', bbox_inches='tight')
    document = '<html><img src="data:image/png;base64,%s"></html>' % base64.b64encode(image_bytes.getvalue()).decode('utf-8')
    debugger.display_html(document, position=2)

def plot_vec(y,length,size=8):
    y = debugger.unwrap(y)
    if y.TypeIsPointerType():
            image_addr = y.GetValueAsUnsigned()
    else:
        image_addr = y.AddressOf().GetValueAsUnsigned()
    data = lldb.process.ReadMemory(image_addr, int(length)*size, lldb.SBError())
    data = np.frombuffer(data, dtype=np.double).flatten()
    x = np.arange(len(data))
    plt.plot(x,data)
    plt.show()
    show()

# if __name__ == "__main__":
#     plt.style.use('ggplot')
#     matplotlib.use( 'tkagg' )
#     plot_vec(np.sin(np.arange(100)*.1))