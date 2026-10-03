"""
Exercice - Simulation de l'acétone (modélisation Python + comparaison DWSIM)

Schéma du procédé (même topologie que dans DWSIM) :

    MSTR-01 --> SPLT-01 --+--> 40 % : reste liquide
                          |
                          +--> 60 % : HEAT-01 --> vapeur saturée

Questions traitées :
    1. température d'ébullition de l'acétone
    2. énergie nécessaire pour vaporiser les 60 %
    3. courbe P = f(T) de 0,1 à 3 atm
"""

import math
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

ZERO_CELSIUS_K = 273.15


# ----------------------------------------------------------------------
# Propriétés physiques de l'acétone
# ----------------------------------------------------------------------
@dataclass(frozen=True)
class Acetone:
    masse_molaire: float = 58.08        # g/mol
    t_critique: float = 508.2           # K

    # Antoine : log10(P_atm) = a - b / (T_K + c)
    a: float = 4.41877
    b: float = 1312.253
    c: float = -32.445

    # Cp liquide (kJ/kg/K) = cp_a + cp_b * T_K
    cp_a: float = 1.337
    cp_b: float = 2.7752e-3

    # ΔHvap (kJ/mol) = dh_a * exp(-beta * Tr) * (1 - Tr)**beta
    dh_a: float = 46.95
    beta: float = 0.2826

    def pression_saturation(self, t_k):
        """Psat (atm) à la température t_k (K)."""
        return 10.0 ** (self.a - self.b / (t_k + self.c))

    def temperature_ebullition(self, p_atm):
        """Température d'ébullition (K) à p_atm : inversion directe d'Antoine
        (corps pur : P = Psat)."""
        return self.b / (self.a - np.log10(p_atm)) - self.c

    def enthalpie_sensible(self, t1_k, t2_k):
        """Intégrale de Cp dT entre t1 et t2 (kJ/kg)."""
        return (self.cp_a * (t2_k - t1_k)
                + 0.5 * self.cp_b * (t2_k ** 2 - t1_k ** 2))

    def enthalpie_vaporisation(self, t_k):
        """Chaleur latente (kJ/kg) à la température t_k."""
        tr = t_k / self.t_critique
        dh_molaire = self.dh_a * math.exp(-self.beta * tr) * (1.0 - tr) ** self.beta
        return dh_molaire * 1000.0 / self.masse_molaire


# ----------------------------------------------------------------------
# Données de l'énoncé
# ----------------------------------------------------------------------
@dataclass
class Alimentation:
    debit_kg_h: float = 30.0            # débit massique de l'énoncé (kg/h)
    t_celsius: float = 10.0
    p_atm: float = 1.0
    part_liquide: float = 0.40
    part_vapeur: float = 0.60


# ----------------------------------------------------------------------
# Calcul du procédé
# ----------------------------------------------------------------------
def simuler(alim: Alimentation, produit: Acetone) -> dict:
    m_liquide = alim.debit_kg_h * alim.part_liquide
    m_vapeur = alim.debit_kg_h * alim.part_vapeur

    t_entree = alim.t_celsius + ZERO_CELSIUS_K
    t_eb = float(produit.temperature_ebullition(alim.p_atm))

    q_sensible = produit.enthalpie_sensible(t_entree, t_eb)
    q_latent = produit.enthalpie_vaporisation(t_eb)
    q_specifique = q_sensible + q_latent

    puissance_kj_h = m_vapeur * q_specifique

    return {
        "m_liquide": m_liquide,
        "m_vapeur": m_vapeur,
        "ecart_bilan": alim.debit_kg_h - (m_liquide + m_vapeur),
        "t_eb_k": t_eb,
        "t_eb_c": t_eb - ZERO_CELSIUS_K,
        "q_sensible": q_sensible,
        "q_latent": q_latent,
        "q_specifique": q_specifique,
        "puissance_kj_h": puissance_kj_h,
        "puissance_kw": puissance_kj_h / 3600.0,
    }


def courbe_saturation(produit: Acetone, p_min=0.1, p_max=3.0, n=25):
    """Retourne (pressions en atm, températures d'ébullition en °C)."""
    pressions = np.linspace(p_min, p_max, n)
    temperatures = produit.temperature_ebullition(pressions) - ZERO_CELSIUS_K
    return pressions, temperatures


# ----------------------------------------------------------------------
# Affichage console
# ----------------------------------------------------------------------
def ligne(libelle, valeur, unite="", decimales=2):
    print(f"  {libelle:<48}{valeur:>12.{decimales}f} {unite}")


def afficher_bilan(alim: Alimentation, r: dict):
    barre = "=" * 72
    print(f"\n{barre}\n  VAPORISATION PARTIELLE DE L'ACÉTONE\n{barre}")

    print("\n  [Alimentation MSTR-01]")
    ligne("Débit massique", alim.debit_kg_h, "kg/h")
    ligne("Température", alim.t_celsius, "°C")
    ligne("Pression", alim.p_atm, "atm")

    print("\n  [Séparateur SPLT-01]")
    ligne("Sortie 1 : 40 % (liquide)", r["m_liquide"], "kg/h")
    ligne("Sortie 2 : 60 % (vers HEAT-01)", r["m_vapeur"], "kg/h")
    ligne("Écart du bilan matière", r["ecart_bilan"], "kg/h", 6)

    print("\n  [Équilibre liquide-vapeur : acétone pure, P = Psat]")
    ligne("Température d'ébullition", r["t_eb_c"], "°C")
    ligne("Température d'ébullition", r["t_eb_k"], "K")

    print("\n  [Réchauffeur HEAT-01]")
    ligne("Chaleur sensible", r["q_sensible"], "kJ/kg")
    ligne("Chaleur latente", r["q_latent"], "kJ/kg")
    ligne("Chaleur totale", r["q_specifique"], "kJ/kg")
    ligne("Puissance pour vaporiser les 60 %", r["puissance_kj_h"], "kJ/h")
    ligne("Puissance pour vaporiser les 60 %", r["puissance_kw"], "kW", 3)
    print(barre)


def afficher_tableau(pressions, temperatures):
    print("\n  n    P (atm)    Teb (°C)")
    print("  " + "-" * 28)
    for i, (p, t) in enumerate(zip(pressions, temperatures), start=1):
        print(f"  {i:<4} {p:7.4f}    {t:8.2f}")


# ----------------------------------------------------------------------
# Graphiques
# ----------------------------------------------------------------------
def graphe_pression_temperature(produit, p_1atm_c, fichier="P_en_fonction_de_T.png"):
    """Graphe demandé : P = f(T) entre 0,1 et 3 atm."""
    p_fin, t_fin = courbe_saturation(produit, n=300)
    p_pts, t_pts = courbe_saturation(produit, n=25)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(t_fin, p_fin, color="#1f5fa8", linewidth=2.2, label="Psat de l'acétone")
    ax.plot(t_pts, p_pts, "s", color="#1f5fa8", markersize=4)
    ax.plot(p_1atm_c, 1.0, "o", color="crimson", markersize=9,
            label=f"Ébullition normale : {p_1atm_c:.1f} °C à 1 atm")
    ax.axhline(1.0, color="grey", linestyle=":", linewidth=1)
    ax.axvline(p_1atm_c, color="grey", linestyle=":", linewidth=1)

    ax.set_xlabel("Température (°C)")
    ax.set_ylabel("Pression (atm)")
    ax.set_title("Acétone : pression de vapeur saturante en fonction de la température")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="upper left")
    fig.tight_layout()
    fig.savefig(fichier, dpi=200)
    return fig


def graphe_ebullition_pression(produit, p_1atm_c, fichier="Teb_en_fonction_de_P.png"):
    """Même courbe, axes inversés : Teb = f(P)."""
    p_fin, t_fin = courbe_saturation(produit, n=300)

    fig, ax = plt.subplots(figsize=(9, 6))
    ax.plot(p_fin, t_fin, color="#2a9d8f", linewidth=2.2, label="Teb de l'acétone")
    ax.plot(1.0, p_1atm_c, "o", color="crimson", markersize=9, label="1 atm")

    ax.set_xlabel("Pression (atm)")
    ax.set_ylabel("Température d'ébullition (°C)")
    ax.set_title("Acétone : température d'ébullition en fonction de la pression")
    ax.grid(True, alpha=0.35)
    ax.legend(loc="lower right")
    fig.tight_layout()
    fig.savefig(fichier, dpi=200)
    return fig


# ----------------------------------------------------------------------
# Programme principal
# ----------------------------------------------------------------------
def main():
    produit = Acetone()
    alimentation = Alimentation()

    resultats = simuler(alimentation, produit)
    afficher_bilan(alimentation, resultats)

    pressions, temperatures = courbe_saturation(produit)
    afficher_tableau(pressions, temperatures)

    graphe_pression_temperature(produit, resultats["t_eb_c"])
    graphe_ebullition_pression(produit, resultats["t_eb_c"])
    print("\n  Graphes enregistrés : P_en_fonction_de_T.png, Teb_en_fonction_de_P.png")
    plt.show()


if __name__ == "__main__":
    main()
