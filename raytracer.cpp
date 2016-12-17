#include <armadillo>
#include <string>
#include <map>
#include <vector>
#include <math.h>

using namespace std;
using namespace arma;

class Material {
public:
	string _nome;
	double _ns, _ni, _d, _sharpness, _illum;
	vec _kd, _ks, _ke, _ka;

	Material() {
		this->_nome = "default";
		this->_ns = 80.0;
		this->_sharpness = 60.0;
		this->_ni = 1.0;
		this->_d = 0.0;
		this->_illum = 1.0;

		this->_kd << 0.5 << 0.5 << 0.5;
		this->_ks << 0.1 << 0.1 << 0.1;
		this->_ke << 0.0 << 0.0 << 0.0;
		this->_ka << 0.0 << 0.0 << 0.0;
	};
}; // Fim do Material

class Luz {
public:
	vec _posicao;
	vec _i;

	Luz(vec posicao) {
		this->_i << 1.0 << 1.0 << 1.0; // Padrão luz branca
		this->_posicao = posicao;
	};
}; // Fim da Luz

class Raio {
public:
	vec _origem, _direcao;
	double _t;

	Raio(vec origem, vec direcao) {
		this->_origem = origem;
		this->_direcao = direcao;
		this->_t = -1;
	};
}; // Fim do Raio

class Face {
public:
	vec A, B, C, _n, BA, CA;
	Material *_mtl;
	double _aux;

	Face(Material *mtl) {
		this->_mtl = mtl;
		this->_n << 0 << 0 << 0;
	};

	virtual double colidir(Raio *raio) {
		return -1;
	};

	Material* material() {
		return this->_mtl;
	};

	virtual vec normal() {
		return this->_n;
	};
}; // Fim da Face

class Esfera : public Face {
public:
	Esfera(Material *mtl, vec centro, double raio) : Face(mtl) {
		A = B = C = _n = centro; // Usar como o centro
		_aux = raio;
	};

	double colidir(Raio *raio) {
		vec distancia = (raio->_origem - this->_n);
		double a = dot(raio->_direcao, raio->_direcao);
		double b = dot(2.0 * raio->_direcao, distancia);
		double c = dot(distancia, distancia) - (this->_aux * this->_aux);
		double delta = (b * b) - 4 * a * c;

		if ((delta > 0)) {
			double raiz = sqrt(delta);
			double dA = (2.0 * a);
			double mB = (-1.0 * b);

			double raiz1 = (mB + raiz) / dA;
			double raiz2 = (mB + raiz) / dA;

			if ((raiz1 > 0) && (raiz1 < raiz2)) {
				return raiz1;
			}

			else if ((raiz2 > 0)) {
				return raiz2;
			}
		}
		return -1;
	};
}; // Fim da Esfera

class Triangulo : public Face {
public:
	Triangulo(Material *mtl, vec v1, vec v2, vec v3) : Face(mtl) {
		this->A = v1;
		this->B = v2;
		this->C = v3;

		this->BA = (this->A - this->B);
		this->CA = (this->A - this->C);
	};

	Triangulo(Material *mtl, vec v1, vec v2, vec v3, vec n) : Triangulo(mtl, v1, v2, v3) {
		this->_n = n;
	};

	vec normal() {
		if ((this->_n.size() <= 0)) {
			this->_n = cross(BA, CA);
		}

		return this->_n;



	}

	// Novo método colidir otimizado, aproveitando operações
	double colidir(Raio *raio) {
		vec OA = (A - raio->_origem);
		mat m(3, 3);

		m.col(0) = this->BA;
		m.col(1) = this->CA;
		m.col(2) = raio->_direcao;
		double dA = det(m);

		if ((dA == 0)) {
			return -1;
		}

		m.col(2) = OA;
		double distancia = (det(m) / dA);

		if ((distancia < 0)) {
			return -1;
		}

		m.col(1) = OA;
		m.col(2) = raio->_direcao;
		double lambda = (det(m) / dA);

		if ((lambda < 0) || (lambda > 1)) {
			return -1;
		}

		m.col(0) = OA;
		m.col(1) = this->CA;
		double beta = (det(m) / dA);

		if ((beta < 0) || (beta > (1 - lambda))) {
			return -1;
		}

		return distancia;
	};

}; // Fim do Triangulo

class Objeto {
public:
	vector<Face*> _faces;
}; // Fim do Objeto

class Cena {
private:
	int _largura, _altura, _total;
	vector<vec> _vetores, _normais;
	map<string, Material*> _materiais;
	vector<Luz*> _luzes;
	vector<Objeto> _objetos;
	vec _camera;
	mat _k;

public:
	Cena() {
		// Objeto padrão caso não haja nenhum "o" no arquivo de entrada
		Objeto obj = Objeto();
		this->_objetos.push_back(obj);

		// Material default
		Material *mtl = new Material();
		this->_materiais[mtl->_nome] = mtl;

		// Padrão caso não tenha configurações customizadas no .obj
		_camera << 0 << 0 << 0;
		_largura = _altura = 100;
	};

	void lerMtl(string entrada) {
		ifstream obj(entrada);
		string linha, identificador;
		Material *mtl = new Material();

		while (getline(obj, linha)) {
			if (!((linha == "") || (linha[0] == '#'))) {
				float a, b, c;
				vec vetor;
				istringstream stream(linha);
				stream >> identificador;

				if ((identificador == "newmtl")) {
					mtl = new Material();
					stream >> mtl->_nome;
					this->_materiais[mtl->_nome] = mtl;
				}
				else if ((identificador == "Ns")) {
					sscanf(linha.c_str(), "%*s %f", &a);
					mtl->_ns = a;
				}
				else if ((identificador == "Ni")) {
					sscanf(linha.c_str(), "%*s %f", &a);
					mtl->_ni = a;
				}
				else if ((identificador == "d")) {
					sscanf(linha.c_str(), "%*s %f", &a);
					mtl->_d = a;
				}
				else if ((identificador == "sharpness")) {
					sscanf(linha.c_str(), "%*s %f", &a);
					mtl->_sharpness = a;
				}
				else if ((identificador == "Ka")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					vetor << a << b << c;
					mtl->_ka = vetor;
				}
				else if ((identificador == "Kd")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					vetor << a << b << c;
					mtl->_kd = vetor;
				}
				else if ((identificador == "Ks")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					vetor << a << b << c;
					mtl->_ks = vetor;
				}
				else if ((identificador == "Ke")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					vetor << a << b << c;
					mtl->_ke = vetor;
				}
			}
		}

		obj.close();
	};

	void lerObj(char* entrada) {
		cout << "Lendo entrada '" << entrada << "'" << endl;

		ifstream obj(entrada);
		string linha, identificador;
		Material *mtl = NULL;
		int faces = 0;

		while (getline(obj, linha)) {
			if (!((linha == "") || (linha[0] == '#'))) {
				float a, b, c;
				vec vetor;

				istringstream stream(linha);
				stream >> identificador;

				if ((identificador == "o")) {
					Objeto objeto = Objeto();
					this->_objetos.push_back(objeto);
				}
				else if ((identificador == "mtllib")) {
					stream >> identificador;
					this->lerMtl(identificador);
				}
				else if ((identificador == "usemtl")) {
					stream >> identificador;
					mtl = this->_materiais[identificador];

					if ((mtl == NULL)) {
						mtl = this->_materiais["default"];
					}
				}
				else if ((identificador == "v")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					vetor << a << b << c;
					this->_vetores.push_back(vetor);
				}
				else if ((identificador == "vn")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					vetor << a << b << c;
					this->_normais.push_back(vetor);
				}
				else if ((identificador == "f")) {
					/*
					* Tenta ler as diversas formas de montar uma face
					* f x y z					> Apenas os vertices
					* f x//a y//a z//a			> Vértices e normais
					* f x/s1/a y/s2/a z/s3/a	> Vértices, superfícies e normais
					*/

					int v1, v2, v3, m, n;
					Triangulo *face = NULL;

					int matches = sscanf(linha.c_str(), "%*s %i %i %i", &v1, &v2, &v3);
					if ((matches == 3)) {
						// Apenas os vértices
						face = new Triangulo(mtl, this->_vetores[v1 - 1], this->_vetores[v2 - 1], this->_vetores[v3 - 1]);
					}
					else {
						matches = sscanf(linha.c_str(), "%*s %i//%i %i//%i %i//%i", &v1, &n, &v2, &n, &v3, &n);
						if ((matches == 6)) {
							// Vértices e normais
							face = new Triangulo(mtl, this->_vetores[v1 - 1], this->_vetores[v2 - 1], this->_vetores[v3 - 1], this->_normais[n - 1]);
						}
						else {
							matches = sscanf(linha.c_str(), "%*s %i/%i/%i %i/%i/%i %i/%i/%i", &v1, &m, &n, &v2, &m, &n, &v3, &m, &n);
							if ((matches == 9)) {
								// Vértices, superfícies(que não utilizamos) e normais
								face = new Triangulo(mtl, this->_vetores[v1 - 1], this->_vetores[v2 - 1], this->_vetores[v3 - 1], this->_normais[n - 1]);
							}
						}
					}

					if ((face != NULL)) {
						++faces;
						Face *f = static_cast<Face*> (face);
						this->_objetos.back()._faces.push_back(f);
						// cout << v1 << ", " << v2 << ", " << v3 << endl;
					}
				}


				/*
				 * Opções customizadas, pra não precisar ficar recompilando sempre...
				 * A partir daqui não existe na documentação do wavefront...
				 */
				else if ((identificador == "light")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					vetor << a << b << c;

					Luz *luz = new Luz(vetor);
					this->_luzes.push_back(luz);
				}
				else if ((identificador == "lightcolor")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					Luz *luz = _luzes.back();
					luz->_i << a << b << c;
				}
				else if ((identificador == "scene")) {
					sscanf(linha.c_str(), "%*s %d %d", &this->_largura, &this->_altura);
				}
				else if ((identificador == "camera")) {
					sscanf(linha.c_str(), "%*s %f %f %f", &a, &b, &c);
					this->_camera << a << b << c;
				}
			}
		}

		obj.close();
		cout << "Pronto! " << this->_objetos.size() << " objeto(s) para renderizar..." << endl;

		// Inicia a matriz K
		this->_k << (this->_largura / 3.0) << 0 << (this->_largura / 2.0) << endr
			<< 0 << -(this->_altura / 3.0) << (this->_altura / 2.0) << endr
			<< 0 << 0 << 1.0;

		this->_k = this->_k.i();
	};

	vec iluminacao(Raio *raio, Face *cFace) {
		// cout << "Que haja luz!!!!" << endl;
		// Todos esses calculos estão no livro, a partir da página 82...

		vec N = normalise(cFace->normal());
		vec kd, ks;
		kd << 0 << 0 << 0; // Lambertian
		ks << 0 << 0 << 0; // Blinn-Phong

		for (int l = 0; l < _luzes.size(); ++l) {
			Luz *luz = _luzes.at(l);
			vec L = normalise(luz->_posicao - raio->_direcao);
			double Lambertian = max(0.0, dot(L, N));

			// Lambertian <= 0: Atrás do objeto
			if ((Lambertian > 0.0)) {
				// raio->posicao - raio->direcao
				// vec H = normalise(_observador + L); // Fórmula do livro
				vec observador = (raio->_origem - raio->_direcao);
				vec H = normalise(observador + L); // Fórmula do livro
				double BlinnPhong = pow(max(0.0, dot(N, H)), cFace->material()->_ns);

				// cout << "Ns: " << cFace->material()->_ns << ", Phong: " << BlinnPhong << endl;

				kd = kd + (cFace->material()->_kd % (luz->_i * Lambertian));
				ks = ks + (cFace->material()->_ks % (luz->_i * BlinnPhong));
			}
		} // Fim das luzes

		vec cor = (cFace->material()->_ka + kd + ks);
		for (int j = 0; j < cor.size(); ++j) {
			cor[j] = (cor[j] >= 1.0) ? 255.0 : (cor[j] <= 0.0) ? 0.0 : (cor[j] * 255.0);
		}

		return cor;
	}

	void renderizar() {
		cout << "Renderizando imagem! " << endl << endl;
		cout << "[ 0% ]\t";

		vec projetivo;
		FILE *img = fopen("imagem.ppm", "w");
		fprintf(img, "P3\n%i %i\n\n255\n", _largura, _altura);
		int x = 0;

		for (int y = 0; y < _altura; ++y) {
			for (x = 0; x < _largura; ++x) {
				projetivo << x << y << 1;
				projetivo = _k * projetivo;

				Raio *raio = new Raio(_camera, projetivo);

				double distancia = -1;
				Face *colisao;

				for (int i = 0; i < _objetos.size(); ++i) {
					Objeto obj = _objetos.at(i);

					for (int j = 0; j < obj._faces.size(); ++j) {
						Face *face = obj._faces.at(j);
						double t = face->colidir(raio);

						if ((t > 0) && ((distancia < 0) || (distancia > t))) {
							distancia = t;
							colisao = face; // Pegar os vetores e o material
						}
					}
				}

				if ((distancia > 0)) {
					vec cor = iluminacao(raio, colisao);
					fprintf(img, "%0.0f %0.0f %0.0f\n", cor[0], cor[1], cor[2]);
				}
				else {
					fprintf(img, "255 255 255\n");
				}
			}

			double total = 1 + ((x * y * 100) / (_altura * _largura));
			if ((total != _total)) {
				_total = total;
				cout << "[ " << _total << "% ]\t";
			}
		}

		fclose(img);
		cout << endl << endl << "Pronto!" << endl;
	};
}; // Fim da Cena

int main(int argc, char **argv) {
	cout << "Iniciando o Ray Tracer..." << endl;

	Cena *cena = new Cena();
	cena->lerObj("entradas/objetos.obj");
	cena->renderizar();

	system("PAUSE");
	return 0;
}