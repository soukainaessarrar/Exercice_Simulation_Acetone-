# Exercice - Simulation de l'acétone

Alimentation : acétone, 10 °C, 1 atm, 30 kg/h. Séparation 40 % / 60 %, les 60 % sont chauffés jusqu'à vaporisation.

## Exécution
```
pip install -r requirements.txt
python exercice_simulation_acetone.py
```

## Sorties
- Bilan matière, température d'ébullition, chaleurs sensible et latente, puissance (kJ/h et kW)
- `P_en_fonction_de_T.png` : P = f(T), de 0,1 à 3 atm
- `Teb_en_fonction_de_P.png` : Teb = f(P)

## DWSIM
Composé : acétone ; modèle NRTL ou Peng-Robinson ; MSTR-01 -> SPLT-01 (0,4 / 0,6) -> HEAT-01 (fraction vapeur = 1).
