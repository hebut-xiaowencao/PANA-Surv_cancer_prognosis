import matplotlib.pyplot as plt
from tensorboard.backend.event_processing.event_accumulator import EventAccumulator

def plot_loss_curve(log_dir):

    event_acc = EventAccumulator(log_dir)
    event_acc.Reload()

    # scalar （cvae_loss、BCE、KLD）
    cvae_loss = event_acc.Scalars('cvae_loss')
    BCE_loss = event_acc.Scalars('BCE')
    KLD_loss = event_acc.Scalars('KLD')

    epochs = [x.step for x in cvae_loss]
    cvae_values = [x.value for x in cvae_loss]
    BCE_values = [x.value for x in BCE_loss]
    KLD_values = [x.value for x in KLD_loss]

    # WCVAE Loss 
    plt.figure(figsize=(8, 5))
    plt.plot(epochs, cvae_values, label='WCVAE Loss', color='b')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('WCVAE Loss Curve')
    plt.legend()
    plt.grid(True)
    plt.show()

    # BCE Loss 
    plt.figure(figsize=(8, 5))
    plt.plot(epochs, BCE_values, label='BCE Loss', color='g')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('BCE Loss Curve')
    plt.legend()
    plt.grid(True)
    plt.show()

    # KLD Loss 
    plt.figure(figsize=(8, 5))
    plt.plot(epochs, KLD_values, label='KLD Loss', color='r')
    plt.xlabel('Epoch')
    plt.ylabel('Loss')
    plt.title('KLD Loss Curve')
    plt.legend()
    plt.grid(True)
    plt.show()

