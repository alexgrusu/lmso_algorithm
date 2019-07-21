#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: Alexandru - George Rusu (2019). All Rights Reserved.
"""
import numpy as np

class cir_buf:
    def __init__(self, cblen):
        self.cblen = cblen
        self.cb = np.zeros([cblen, 1])
        
    def cb_push(self, x):
        self.cb = np.vstack((x, self.cb[0 : self.cblen - 1]))
    
    def dot_product(self):
        return ((self.cb).T).dot(self.cb)


class fir(cir_buf):
    def __init__(self, ford):
        self.ord = ford
        self.coeffs = np.zeros([ford, 1])
    
    def ffir(self, x):
        return ((x.cb).T).dot(self.coeffs)


class lmso(fir, cir_buf):
  def __init__(self, ford = 512, cst = 2.165e-6, lam = 0.95, sgmv = 0, sgmw = 0, lniproc = 2, alpha = 1.0, delta = 1e-6):
    id_vec = np.ones([ford, 1])
    id_matrix = np.diagflat(id_vec)
    
    self.h = fir(ford)
    self.c = np.sqrt(cst) * id_vec
    self.ccorr_mat = cst * id_matrix
    self.gamma_mat = self.ccorr_mat + sgmw * id_matrix
    
    self.x = cir_buf(ford)
    self.xcorr_mat = np.zeros([ford, ford])
    
    self.msd = cst
    self.tmp = 0
    self.sx = 0
    self.sv = sgmv
    self.sw = sgmw
    self.sw_idmat = sgmw * id_matrix
    
    self.ln_init_proc = lniproc * ford
    self.alpha = alpha
    self.delta = delta
    self.lam = lam
   
  def lmso_w(self, echo):
    if(self.ln_init_proc == 0):
      self.msd = np.trace(self.ccorr_mat) + self.h.ord * self.sw
      self.ln_init_proc = -1
    
    err = echo - self.h.ffir(self.x)
    if(self.ln_init_proc > 0):
       self.sx = self.x.dot_product()
       sz = self.alpha / (self.sx + self.delta)
       u = sz * err * self.x.cb
       self.c += u
       self.ccorr_mat = self.lam * self.ccorr_mat + (1 - self.lam) * (self.c).dot((self.c).T)
       self.ln_init_proc -= 1
    else:
       self.sx = self.x.dot_product() / self.h.ord
       sz = self.msd / (self.msd * self.sx * (self.h.ord + 2) + self.h.ord * self.sv)
       u = sz * err * self.x.cb
       self.tmp = self.msd * (1 - sz * self.sx)
       self.msd = self.tmp + self.h.ord * self.sw
    
    self.h.coeffs += u
    return err

  def lmso_g(self, echo):
    self.xcorr_mat = self.lam * self.xcorr_mat + (1 - self.lam) * (self.x.cb).dot((self.x.cb).T)
    err = echo - self.h.ffir(self.x)
    
    if(self.ln_init_proc > 0):
       self.sx = self.x.dot_product()
       sz = self.alpha / (self.sx + self.delta)
       u = sz * err * self.x.cb
       self.c += u
       self.ccorr_mat = self.lam * self.ccorr_mat + (1 - self.lam) * (self.c).dot((self.c).T)
       self.ln_init_proc -= 1
    else:
       self.sx = self.x.dot_product() / self.h.ord
       gr_prod = (self.gamma_mat).dot(self.xcorr_mat)
       rg_prod = (self.xcorr_mat).dot(self.gamma_mat)
       tp = np.trace(gr_prod)
       sz = tp / (self.h.ord * self.sx * tp + 2 * np.trace(gr_prod.dot(self.xcorr_mat)) \
                  + self.h.ord * self.sx * self.sv)
       u = sz * err * self.x.cb
       self.ccorr_mat = self.gamma_mat - sz * (gr_prod + rg_prod) \
                      + sz**2 * (2 * rg_prod.dot(self.xcorr_mat) + self.xcorr_mat * tp) \
                      + sz**2 * self.sv * self.xcorr_mat
    
    self.h.coeffs += u
    self.gamma_mat = self.ccorr_mat + self.sw_idmat
    return err

