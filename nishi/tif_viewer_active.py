import cv2
from IPython.display import Image, display
from ipywidgets import widgets


def imshow(img):
    """画像を Notebook 上に表示する。
    """
    encoded = cv2.imencode(".png", img)[1]
    display(Image(encoded, width=400))


def process(img, alpha, beta):
    """明るさ、コントラストを調整し、結果を表示する。
    """
    dst = adjust(img, alpha, beta)
    imshow(dst)


param_widgets = {}

# パラメータ「ゲイン」を設定するスライダー
param_widgets["alpha"] = widgets.FloatSlider(
    min=0.0, max=3.0, step=0.1, value=1.0, description="alpha: "
)

# パラメータ「バイアス」を設定するスライダー
param_widgets["beta"] = widgets.FloatSlider(
    min=-100.0, max=100.0, step=10.0, value=0.0, description="beta: "
)

for x in param_widgets.values():
    x.layout.width = "400px"

# 画像を読み込む。
img = cv2.imread("File74.tif")

# ウィジェットを表示する。
widgets.interactive(process, img=widgets.fixed(img), **param_widgets)