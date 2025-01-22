from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtWidgets import QSpinBox, QGraphicsDropShadowEffect
from PyQt5.QtGui import QColor, QPalette

__all__ = ['changeSkinDark', 'changeSkinLight']

# COLOR THEMES
@QtCore.pyqtSlot()
def changeSkinDark(self):
    # PALETTE
    dark_palette = QPalette()
    
    
    # Active
    dark_palette.setColor(QPalette.Window, QColor(35, 35, 40))
    dark_palette.setColor(QPalette.WindowText, QColor(200,200,200))
    dark_palette.setColor(QPalette.Base, QColor(60, 60, 65))
    dark_palette.setColor(QPalette.AlternateBase, QColor(35, 35, 40))
    #dark_palette.setColor(QPalette.ToolTipBase, QColor(255,255,255))
    #dark_palette.setColor(QPalette.ToolTipText, QColor(255,255,255))
    dark_palette.setColor(QPalette.Text, QColor(200,200,200))
    dark_palette.setColor(QPalette.Button, QColor(35, 35, 40))
    dark_palette.setColor(QPalette.ButtonText, QColor(200,200,200))
    #dark_palette.setColor(QPalette.BrightText, QColor(255,0,0))
    #dark_palette.setColor(QPalette.Link, QColor(42, 130, 218))
    dark_palette.setColor(QPalette.Highlight, QColor(42, 130, 218))
    
    # Disabled
    dark_palette.setColor(QPalette.Disabled, QPalette.ButtonText, QColor(100,100,100))

    QtGui.QGuiApplication.setPalette(dark_palette)
    
    if hasattr(self,"INACTIVE_SOLVENT_TAB_2"):
        self.INACTIVE_SOLVENT_TAB_2.setStyleSheet("QLabel { background-color: rgb(47, 47, 52)}")
    
    # Radio Buttons Customization
    size = self.style().pixelMetric(QtWidgets.QStyle.PM_ExclusiveIndicatorWidth)
    border = 1
    for rb in [self.extremes_radioButton,self.centerRange_radioButton]:
        rb.setStyleSheet('''
            QRadioButton::indicator {{
                border: {border}px solid rgb(26, 26, 30); 
                height: {size}px;
                width: {size}px;
                border-radius: {radius}px;
            }}
            QRadioButton::indicator:checked {{
                background: qradialgradient(
                    cx:.5, cy:.5, radius: {innerRatio},
                    fx:.5, fy:.5,
                    stop:0 {checkColor}, 
                    stop:0.29 {checkColor},
                    stop:0.3 rgb(60, 60, 65),
                    stop:1 rgb(60, 60, 65)
                    );
            }}
            QRadioButton::indicator:focus {{
                outline: none;
            }}
        '''.format(
            size=size - border * 2, 
            border=border, 
            radius=size // 2, 
            innerRatio=1 - (border * 2 + 1) / size, 
            checkColor='rgb(42, 130, 218)'
        ))
    
    # Measurement Tab Charts
    self.rms_text.txt.get_children()[0].set_color("white")
    self.rms_text.patch.set_facecolor((60/255,60/255,65/255))
    self.rms_text.patch.set_edgecolor((200/255,200/255,200/255,1))

    for chart in self.charts.values():
        chart.fig.patch.set_facecolor((35/255,35/255,40/255,1))
        chart.axes.set_facecolor((35/255,35/255,40/255,1))
        
        for spine in chart.axes.spines.values():
            spine.set_color((200/255,200/255,200/255,1))
        
        chart.axes.set_xlabel(chart.axes.get_xlabel(), fontdict={'color': (200/255,200/255,200/255,1)})
        chart.axes.set_ylabel(chart.axes.get_ylabel(), fontdict={'color': (200/255,200/255,200/255,1)})
        chart.axes.tick_params(axis='both',which='both',colors=(200/255,200/255,200/255,1))
        chart.axes.grid(which="both",color=(60/255,60/255,65/255,1))
        chart.draw_idle()
    
    # Fitting Tab Charts
    for chart_types in self.fitting_charts.values():
        for chart in chart_types.values():
            chart.fig.patch.set_facecolor((35/255,35/255,40/255,1))
            chart.axes.set_facecolor((35/255,35/255,40/255,1))
            
            for spine in chart.axes.spines.values():
                spine.set_color((200/255,200/255,200/255,1))
            
            chart.axes.set_title(chart.axes.get_title(), fontdict={'color': (200/255,200/255,200/255,1)})
            chart.axes.set_xlabel(chart.axes.get_xlabel(), fontdict={'color': (200/255,200/255,200/255,1)})
            chart.axes.set_ylabel(chart.axes.get_ylabel(), fontdict={'color': (200/255,200/255,200/255,1)})
            chart.axes.tick_params(axis='both',which='both',colors=(200/255,200/255,200/255,1))
            chart.axes.grid(which="both",color=(60/255,60/255,65/255,1))
            chart.draw_idle()
    
    # STYLESHEET
    stylesheet = """ 
        QTabBar::tab {height: 1em; margin: 0px; padding: 4px; padding-left: 1em; padding-right: 1em; /* height: 1em is expected to recover original setting */
            border-top: 1px solid rgba(60, 60, 65, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
            border-left: 1px solid rgba(60, 60, 65, 1);
            border-right: 1px solid rgba(60, 60, 65, 1);}
        QTabBar::tab:!selected {margin-top: 3px;}
        QTabBar::tab:selected {height: 15 px; border: 1px solid rgba(42, 130, 218, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
            background: rgba(42, 130, 218, 1); color: white} /* height: 15 px is expected to cover for slight change */
                                                                /* in the height when only one tab is present in the tab bar. */
                                                                /* It is prone to font size changes though. */
        """

    self.setStyleSheet(stylesheet)

@QtCore.pyqtSlot()
def changeSkinLight(self):
    # PALETTE
    QtGui.QGuiApplication.setPalette(self.default_palette)
    
    if hasattr(self,"INACTIVE_SOLVENT_TAB_2"):
        self.INACTIVE_SOLVENT_TAB_2.setStyleSheet("QLabel { background-color: rgb(252, 252, 252)}")
    
    # Radio Buttons Customization
    size = self.style().pixelMetric(QtWidgets.QStyle.PM_ExclusiveIndicatorWidth)
    border = 1
    for rb in [self.extremes_radioButton,self.centerRange_radioButton]:
        rb.setStyleSheet('''
            QRadioButton::indicator {{
                border: {border}px solid rgb(26, 26, 30); 
                height: {size}px;
                width: {size}px;
                border-radius: {radius}px;
            }}
            QRadioButton::indicator:checked {{
                background: qradialgradient(
                    cx:.5, cy:.5, radius: {innerRatio},
                    fx:.5, fy:.5,
                    stop:0 {checkColor}, 
                    stop:0.29 {checkColor},
                    stop:0.3 rgb(255,255,255),
                    stop:1 rgb(255,255,255)
                    );
            }}
            QRadioButton::indicator:focus {{
                outline: none;
            }}
        '''.format(
            size=size - border * 2, 
            border=border, 
            radius=size // 2, 
            innerRatio=1 - (border * 2 + 1) / size, 
            checkColor='rgb(42, 130, 218)'
        ))
    
    # Measurement Tab Charts
    self.rms_text.txt.get_children()[0].set_color("black")
    self.rms_text.patch.set_facecolor("white")
    self.rms_text.patch.set_edgecolor("black")
    
    for chart in self.charts.values():
        chart.fig.patch.set_facecolor((255/255,255/255,255/255,1))
        chart.axes.set_facecolor((255/255,255/255,255/255,1))

        for spine in chart.axes.spines.values():
            spine.set_color((0,0,0,1))

        chart.axes.set_xlabel(chart.axes.get_xlabel(), fontdict={'color': (0,0,0,1)})
        chart.axes.set_ylabel(chart.axes.get_ylabel(), fontdict={'color': (0,0,0,1)})
        chart.axes.tick_params(axis='both',which='both',colors=(0,0,0,1))
        chart.axes.grid(which="both",color='#f9f9f9')
        chart.draw_idle()
    
    # Fitting Tab Charts
    for chart_types in self.fitting_charts.values():
        for chart in chart_types.values():
            chart.fig.patch.set_facecolor((255/255,255/255,255/255,1))
            chart.axes.set_facecolor((255/255,255/255,255/255,1))
            
            for spine in chart.axes.spines.values():
                spine.set_color((0,0,0,1))
            
            chart.axes.set_title(chart.axes.get_title(), fontdict={'color': (0,0,0,1)})
            chart.axes.set_xlabel(chart.axes.get_xlabel(), fontdict={'color': (0,0,0,1)})
            chart.axes.set_ylabel(chart.axes.get_ylabel(), fontdict={'color': (0,0,0,1)})
            chart.axes.tick_params(axis='both',which='both',colors=(0,0,0,1))
            chart.axes.grid(which="both",color='#f9f9f9')
            chart.draw_idle()
    
    # STYLESHEET
    stylesheet = """ 
        QTabBar::tab {height: 1em; margin: 0px; padding: 4px; padding-left: 1em; padding-right: 1em; /* height: 1em is expected to recover original setting */
            border-top: 1px solid rgba(200, 200, 205, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
            border-left: 1px solid rgba(200, 200, 205, 1);
            border-right: 1px solid rgba(200, 200, 205, 1);
            }
        QTabBar::tab:!selected {margin-top: 3px;}
        QTabBar::tab:selected {height: 15px; border: 1px solid rgba(42, 130, 218, 1); border-top-left-radius: 3px; border-top-right-radius: 3px;
            background: rgba(42, 130, 218, 1); color: white} /* height: 15 px is expected to cover for slight change */
                                                                /* in the height when only one tab is present in the tab bar. */
                                                                /* It is prone to font size changes though. */
        """

    self.setStyleSheet(stylesheet)

# OTHER STYLES
@QtCore.pyqtSlot()
def changeReadOnlyStyle(self):
    # STYLESHEET
    stylesheet = """ 
        QDoubleSpinBox[readOnly="false"] { font-style: normal; }
        QDoubleSpinBox[readOnly="true"] { font-style: italic; }
        QComboBox[enabled="false"] { font-style: italic; }
        QComboBox[enabled="true"] { font-style: normal; }
        """
    
    current_stylesheet = self.styleSheet()
    
    self.setStyleSheet(current_stylesheet + stylesheet)

@QtCore.pyqtSlot()
def improperValueStyle(self, widget: QSpinBox):
    current_value = widget.value()
    min_value = widget.minimum()
    max_value = widget.maximum()
    
    color = ""
    
    if current_value == min_value or current_value == max_value:
        color = "red"
    
    # Update the spinbox stylesheet
    if color:
        widget.setStyleSheet(f"""QDoubleSpinBox
                             {{ border: 1px solid {color};
                                border-radius: 2px;
                                margin-top: 2px;
                                margin-bottom: 2px;
                                margin-left: 0px;
                                margin-right: 0px;
                             }}""")
        effect = QGraphicsDropShadowEffect()
        effect.setBlurRadius(20)
        effect.setColor(QColor(255, 0, 0, 128)) # Glow color: red
        effect.setOffset(0, 0) # Apply the effect to the widget
        widget.setGraphicsEffect(effect)
    
    else:
        widget.setStyleSheet("")
        widget.setGraphicsEffect(None)