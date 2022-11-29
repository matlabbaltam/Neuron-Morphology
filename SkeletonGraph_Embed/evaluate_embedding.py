from sklearn import preprocessing
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from sklearn.metrics import accuracy_score
from sklearn.model_selection import GridSearchCV, KFold, StratifiedKFold
from sklearn.model_selection import cross_val_score
from sklearn.neural_network import MLPClassifier
from sklearn.svm import SVC, LinearSVC
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import torch
import torch.nn as nn

def draw_plot(datadir, DS, embeddings, fname, max_nodes=None):
    return
    import seaborn as sns
    graphs = read_graphfile(datadir, DS, max_nodes=max_nodes)
    labels = [graph.graph['label'] for graph in graphs]

    labels = preprocessing.LabelEncoder().fit_transform(labels)
    x, y = np.array(embeddings), np.array(labels)
    print('fitting TSNE ...')
    x = TSNE(n_components=2).fit_transform(x)

    plt.close()
    df = pd.DataFrame(columns=['x0', 'x1', 'Y'])

    df['x0'], df['x1'], df['Y'] = x[:,0], x[:,1], y
    sns.pairplot(x_vars=['x0'], y_vars=['x1'], data=df, hue="Y", size=5)
    plt.legend()
    plt.savefig(fname)

class LogReg(nn.Module):
    def __init__(self, ft_in, nb_classes):
        super(LogReg, self).__init__()
        self.fc = nn.Linear(ft_in, nb_classes)

        for m in self.modules():
            self.weights_init(m)

    def weights_init(self, m):
        if isinstance(m, nn.Linear):
            torch.nn.init.xavier_uniform_(m.weight.data)
            if m.bias is not None:
                m.bias.data.fill_(0.0)

    def forward(self, seq):
        ret = self.fc(seq)
        return ret

def logistic_classify(x, y):
    nb_classes = np.unique(y).shape[0]
    xent = nn.CrossEntropyLoss()
    hid_units = x.shape[1]

    accs = []
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    for train_index, test_index in kf.split(x, y):
        train_embs, test_embs = x[train_index], x[test_index]
        train_lbls, test_lbls= y[train_index], y[test_index]

        train_embs, train_lbls = torch.from_numpy(train_embs).cuda(), torch.from_numpy(train_lbls).cuda()
        test_embs, test_lbls= torch.from_numpy(test_embs).cuda(), torch.from_numpy(test_lbls).cuda()


        log = LogReg(hid_units, nb_classes)
        log.cuda()
        opt = torch.optim.Adam(log.parameters(), lr=0.01, weight_decay=0.0)

        best_val = 0
        test_acc = None
        for it in range(100):
            log.train()
            opt.zero_grad()

            logits = log(train_embs)
            loss = xent(logits, train_lbls)

            loss.backward()
            opt.step()

        logits = log(test_embs)
        preds = torch.argmax(logits, dim=1)
        acc = torch.sum(preds == test_lbls).float() / test_lbls.shape[0]
        accs.append(acc.item())
    return np.mean(accs)

def svc_classify(x, y, search):
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    accuracies = []
    for train_index, test_index in kf.split(x, y):

        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # x_train, x_val, y_train, y_val = train_test_split(x_train, y_train, test_size=0.1)
        if search:
            params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
            classifier = GridSearchCV(SVC(), params, cv=5, scoring='accuracy', verbose=0)
        else:
            classifier = SVC(C=10)
        classifier.fit(x_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(x_test)))
    return np.mean(accuracies)

def randomforest_classify(x, y, search):
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    accuracies = []
    for train_index, test_index in kf.split(x, y):

        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        if search:
            params = {'max_depth':[3,5,7,9],'random_state':[0]}
            classifier = GridSearchCV(RandomForestClassifier(), params, cv=5, scoring='accuracy', verbose=0)
        else:
            classifier = RandomForestClassifier()
        classifier.fit(x_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(x_test)))
    ret = np.mean(accuracies)
    return ret

def linearsvc_classify(x, y, search):
    kf = StratifiedKFold(n_splits=10, shuffle=True, random_state=None)
    accuracies = []
    for train_index, test_index in kf.split(x, y):

        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        if search:
            params = {'C':[0.001, 0.01,0.1,1,10,100,1000]}
            classifier = GridSearchCV(LinearSVC(), params, cv=5, scoring='accuracy', verbose=0)
        else:
            classifier = LinearSVC(C=10)
        classifier.fit(x_train, y_train)
        accuracies.append(accuracy_score(y_test, classifier.predict(x_test)))
    return np.mean(accuracies)


def knncluster(x,y):
    kmeans = KMeans(n_clusters=100, init='k-means++',random_state=0).fit(x)
    cluster_labels = kmeans.labels_
    cluster2label = {}
    for cluster_id in cluster_labels:
        cell_indx = np.where(cluster_labels==cluster_id)[0]
        cell_indx = cell_indx.astype(int)
        #cell_labels = labels[cell_indx]
        tmp_labels = []
        for one in cell_indx:
            tmp_labels.append(y[one])
        label = max(set(tmp_labels), key = tmp_labels.count)
        cluster2label[cluster_id] = label    
    prediction_cluster = []
    for cluster_id in cluster_labels:
        prediction_cluster.append(cluster2label[cluster_id])
    false = 0
    for i in range(len(prediction_cluster)):
        #prediction_str.append(all_cell_type[prediction[i]])
        if prediction_cluster[i]!=y[i]:
            false=false+1
    #print (1-false/len(prediction_cluster))
    return (1-false/len(prediction_cluster))



def evaluate_embedding(embeddings, labels, search=True):
    labels = preprocessing.LabelEncoder().fit_transform(labels)
    x, y = np.array(embeddings), np.array(labels)
    # print(x.shape, y.shape)

    logreg_accuracies = logistic_classify(x, y) 
    # print(logreg_accuracies)
    print('LogReg', np.mean(logreg_accuracies))

    svc_accuracies = svc_classify(x,y, search)
    # print(svc_accuracies)
    print('svc', np.mean(svc_accuracies))

    linearsvc_accuracies = linearsvc_classify(x, y, search) 
    # print(linearsvc_accuracies)
    print('LinearSvc', np.mean(linearsvc_accuracies))

    randomforest_accuracies = randomforest_classify(x, y, search)
    # print(randomforest_accuracies)
    print('randomforest', np.mean(randomforest_accuracies))

    clusteraccuracy = knncluster(x,y)
    print('KNNclustering',np.mean(clusteraccuracy))

    return np.mean(logreg_accuracies), np.mean(svc_accuracies), np.mean(linearsvc_accuracies), np.mean(randomforest_accuracies),np.mean(clusteraccuracy)

if __name__ == '__main__':
    evaluate_embedding('./data', 'ENZYMES', np.load('tmp/emb.npy'))
