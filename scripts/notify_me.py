# -*- coding: utf-8 -*-
"""
Created on Thu Mar 15 17:13:00 2018

Code adapted from: https://stackoverflow.com/a/16509278

@author: Calil
"""

import smtplib
import os.path as op
from email.mime.multipart import MIMEMultipart
from email.mime.base import MIMEBase
from email.mime.text import MIMEText
from email.utils import COMMASPACE, formatdate
from email import encoders

"""
The following function was adapted from: https://stackoverflow.com/a/16509278
"""
def send_mail(send_from, send_to, subject, message, files=[],
              server="localhost", port=587, username='', password='',
              use_tls=True):
    """Compose and send email with provided info and attachments.

    Args:
        send_from (str): from name
        send_to (str): to name
        subject (str): message title
        message (str): message body
        files (list[str]): list of file paths to be attached to email
        server (str): mail server host name
        port (int): port number
        username (str): server auth username
        password (str): server auth password
        use_tls (bool): use TLS mode
    """
    msg = MIMEMultipart()
    msg['From'] = send_from
    msg['To'] = COMMASPACE.join(send_to)
    msg['Date'] = formatdate(localtime=True)
    msg['Subject'] = subject

    msg.attach(MIMEText(message))

    for path in files:
        part = MIMEBase('application', "octet-stream")
        with open(path, 'rb') as file:
            part.set_payload(file.read())
        encoders.encode_base64(part)
        part.add_header('Content-Disposition',
                        'attachment; filename="{}"'.format(op.basename(path)))
        msg.attach(part)

    smtp = smtplib.SMTP(server, port)
    if use_tls:
        smtp.starttls()
    smtp.login(username, password)
    smtp.sendmail(send_from, send_to, msg.as_string())
    smtp.quit()
    
if __name__ == '__main__':
    file = open('bot_email')
    sender = file.readline()
    pwrd = file.readline()
    receiver = file.readline()
    subject = "Your simulation is done :D"
    msg = "The simulation results are attached to this email \o/"
    fls = ["simulation_results.tar"]
    svr = 'smtp.gmail.com'
    user = sender
    file.close()

    send_mail(sender,receiver,subject,msg,files = fls,server = svr, username=user, password=pwrd)